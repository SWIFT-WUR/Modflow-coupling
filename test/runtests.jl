using Test
using Dates
using Statistics
using NCDatasets
using IndusModflow
import IndusModflow as IM
import BasicModelInterface as BMI

# Static inputs: INDUS_MODFLOW_DATA, or `data/` next to this repository
const DATA_DIR = get(ENV, "INDUS_MODFLOW_DATA", normpath(@__DIR__, "..", "..", "data"))

function test_config(dir; start = Date(1979, 1, 1), stop = Date(1979, 3, 31), kwargs...)
    Config(; data_dir = DATA_DIR, work_dir = joinpath(dir, "mf6"), output_dir = joinpath(dir, "output"),
           start_date = start, end_date = stop, kwargs...)
end

"Synthetic VIC forcing: uniform recharge on land and a fraction of bankfull discharge."
function force!(m, frac; recharge = 0.5e-3)
    land = m.inputs.landmask
    BMI.set_value(m, "groundwater_recharge", vec(ifelse.(land, recharge, NaN)))
    BMI.set_value(m, "river_discharge", vec(ifelse.(land, frac .* m.inputs.qbank, NaN)))
end

"Budget rates (in, out) of a package from the MF6 listing file, per stress period."
function listing_rates(path, package)
    blocks = split(read(path, String), "VOLUME BUDGET FOR ENTIRE MODEL")[2:end]
    map(blocks) do b
        rates = [parse(Float64, mm[1]) for mm in eachmatch(Regex("$package\\s*=\\s*[-0-9.E+]+\\s+$package\\s*=\\s*([-0-9.E+]+)"), b)]
        (rates[1], rates[2])
    end
end

const FORCING = (0.3, 0.5, -0.1)   # March has negative discharge

@testset "IndusModflow" begin
    tmp = mktempdir()
    cfg = test_config(tmp)
    s = read_static_inputs(cfg)
    p = build_parameters(s, cfg)

    @testset "Case 0: static audit" begin
        for (name, ok) in audit_static_inputs(s, p)
            @test ok
        end
        @test 2 * length(p.chd_cells) == 1410          # reference_cases.md regression count
        @test length(p.river.cells) == 15165
        @test all(p.chd_head .== s.steady_head[p.chd_cells, :])
    end

    @testset "parameter formulas" begin
        @test IM.permeability_to_conductivity(-12.0) ≈ 1000 * 1e-12 * (9.81 / 0.001) * 86400
        r = p.river
        i = findfirst(r.is_river)
        @test IM.river_stage(r, i, -1.0) == r.bottom[i]            # no river infiltration
        @test IM.river_stage(r, i, 0.0) == r.bottom_bankfull[i]    # depth < 0.01 m -> 0
        q = 1000.0
        depth = (0.045q / (r.width[i] * sqrt(r.slope[i])))^(3 / 5)
        @test IM.river_stage(r, i, q) ≈ r.bottom_bankfull[i] + depth
        j = findfirst(!, r.is_river)
        @test IM.river_stage(r, j, 50.0) == r.head_ini[j]
        c = p.evt_cells[1]
        @test p.evt_rate[1] ≈ s.cprate[c] * s.cell_area[c] / (s.delr * s.delc)
    end

    @testset "Case 1/2: continuous run, BMI exchange" begin
        m = IndusGroundwater(cfg)
        @test BMI.get_current_time(m) == 0.0
        @test BMI.get_end_time(m) == 90.0
        @test BMI.get_time_step(m) == 31.0
        @test sort(BMI.get_input_var_names(m)) == ["groundwater_abstraction", "groundwater_recharge", "river_discharge"]
        riv, evt, rch = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
        heads = Array{Float64,3}[]
        for frac in FORCING
            force!(m, frac)
            BMI.update(m)
            push!(riv, package_flows(m, "RIV")); push!(evt, package_flows(m, "EVT")); push!(rch, package_flows(m, "RCH"))
            push!(heads, copy(m.head))
            # MODFLOW -> VIC conversions (coupling_contract.md)
            bf = reshape(BMI.get_value(m, "baseflow_to_routing", zeros(length(m.baseflow))), size(m.baseflow))
            @test all(isnan, bf[.!s.landmask])
            @test all(bf[c] ≈ -riv[end][k] / 86400 for (k, c) in enumerate(p.river.cells))
            days = daysinmonth(m.periods[m.completed])
            @test all(m.capillary[c] ≈ -evt[end][k] * days / s.cell_area[c] * 1000 for (k, c) in enumerate(p.evt_cells))
            @test all(>=(0), m.capillary)
        end
        @test BMI.get_current_time(m) == 90.0
        # RCH volume uses the true cell area, not the uniform DIS area. MF6
        # reports no boundary flow at constant-head cells (idomain == 2).
        land_variable_head = s.landmask .& (s.idomain .== 1)
        @test sum(rch[1]) ≈ 0.5e-3 * sum(s.cell_area[land_variable_head]) rtol = 1e-10
        @test all(<=(0), riv[3])                           # negative discharge: no river infiltration
        @test all(isnan, heads[1][:, :, 1][s.idomain .== 0])
        @test all(heads[3][c, l] ≈ p.chd_head[k, l] for (k, c) in enumerate(p.chd_cells), l in 1:2)
        BMI.finalize(m)

        # BMI SIMVALS equal the MF6 budget in the listing file
        lst = joinpath(cfg.work_dir, "indus.lst")
        for (pkg, flows) in (("RIV", riv), ("EVT", evt), ("RCH", rch))
            for (k, (qin, qout)) in enumerate(listing_rates(lst, pkg))
                @test sum(max.(flows[k], 0)) ≈ qin rtol = 1e-8 atol = 1e-3
                @test -sum(min.(flows[k], 0)) ≈ qout rtol = 1e-8 atol = 1e-3
            end
        end

        # Monthly output file
        NCDataset(joinpath(cfg.output_dir, "indus_modflow_monthly.nc")) do ds
            @test length(ds["time"]) == 3
            @test isequal(coalesce.(ds["groundwater_head"][:, :, :, 3], NaN), heads[3])
        end

        @testset "restart after January equals continuous run" begin
            restart = test_config(joinpath(tmp, "restart"); start = Date(1979, 2, 1),
                                  initial_head_file = joinpath(cfg.output_dir, "head_19790131.nc"))
            m2 = IndusGroundwater(restart)
            for (k, frac) in enumerate(FORCING[2:3])
                force!(m2, frac)
                BMI.update(m2)
                active = s.idomain .> 0
                d = abs.(m2.head .- heads[k + 1])[cat(active, active; dims = 3)]
                @test maximum(d) < 1e-6
                @test package_flows(m2, "RIV") ≈ riv[k + 1] rtol = 1e-6
            end
            BMI.finalize(m2)
        end
    end

    @testset "forcing file and BMI metadata" begin
        include(joinpath(@__DIR__, "..", "scripts", "make_synthetic_forcing.jl"))
        forcing = make_synthetic_forcing(joinpath(tmp, "forcing", "vic.nc"), 1979:1979)
        cfgf = test_config(joinpath(tmp, "forcing"); stop = Date(1979, 2, 28), forcing_file = forcing)
        m = IndusGroundwater(cfgf)
        @test BMI.get_time_units(m) == "days since 1979-01-01"
        @test BMI.get_grid_shape(m, 0) == [180, 204]
        @test BMI.get_grid_shape(m, 1, zeros(Int, 3)) == [2, 180, 204]
        @test BMI.get_grid_origin(m, 0, zeros(2)) == [s.lat[1], s.lon[1]]
        c = p.river.cells[findfirst(p.river.is_river)]
        vic = NCDataset(forcing) do ds   # VIC layout: latitude north-first
            (ds["OUT_GWRECHARGE"][c[1], end - c[2] + 1, 1], ds["OUT_DISCHARGE"][c[1], end - c[2] + 1, 1])
        end
        BMI.update(m)                    # January from the file
        @test m.recharge[c] ≈ vic[1] / 31 / 1000
        @test m.discharge[c] ≈ vic[2]
        BMI.set_value(m, "groundwater_recharge", fill(1e-3, length(m.recharge)))
        BMI.update(m)                    # February: recharge set through BMI wins, discharge from file
        @test all(==(1e-3), m.recharge)
        @test m.discharge[c] != vic[2]
        layer2 = BMI.get_value(m, "groundwater_head_layer_2", zeros(length(m.recharge)))
        @test isequal(layer2, vec(m.head[:, :, 2]))
        depth = BMI.get_value(m, "groundwater_depth", zeros(length(m.recharge)))
        land = findfirst(vec(s.landmask .& (s.idomain .> 0)))
        @test depth[land] ≈ vec(p.top)[land] - vec(m.head[:, :, 1])[land]
        # A new model finalizes the active one instead of clashing in libmf6
        m2 = IndusGroundwater(cfgf)
        @test m.mf === nothing
        BMI.update(m2)
        BMI.finalize(m2)
    end

    @testset "input validation" begin
        cfgw = test_config(joinpath(tmp, "wel"); stop = Date(1979, 1, 31), pumping = true)
        m = IndusGroundwater(cfgw)
        bad = zeros(size(m.abstraction)); bad[findfirst(==(0), s.idomain)] = 1.0
        @test_throws ErrorException BMI.set_value(m, "groundwater_abstraction", vec(bad))
        @test_throws ErrorException BMI.set_value(m, "groundwater_abstraction", fill(-1.0, length(bad)))
        pump = zeros(size(m.abstraction)); cell = p.rch_cells[findfirst(c -> s.idomain[c] == 1, p.rch_cells)]
        pump[cell] = 1e4
        force!(m, 0.3)
        BMI.set_value(m, "groundwater_abstraction", vec(pump))
        BMI.update(m)
        wel = package_flows(m, "WEL")
        @test sum(wel) ≈ -1e4                     # WEL rate is the negative abstraction
        BMI.finalize(m)
        @test_throws ErrorException align_to_grid(zeros(204, 180), s.lon .+ 0.01, s.lat, s)
        shuffled = reverse(s.dem; dims = 2)
        @test isequal(align_to_grid(shuffled, s.lon, reverse(s.lat), s), s.dem)
    end
end
