using Test
using Dates
using Statistics
using NCDatasets
using ModflowCoupling
import ModflowCoupling as MC
import BasicModelInterface as BMI

# The tests run the Indus example model, with monthly periods from the global
# steady state (the Python reference cases). Static inputs: INDUS_DATA, or
# `data/` next to this repository (the data_dir of examples/indus/indus.toml).
const INDUS_DIR = normpath(@__DIR__, "..", "examples", "indus")
const INDUS = load_config(joinpath(INDUS_DIR, "indus.toml"))
const DATA_DIR = get(ENV, "INDUS_DATA", INDUS.data_dir)

function test_config(dir; start = Date(1979, 1, 1), stop = Date(1979, 3, 31), kwargs...)
    Config(INDUS; data_dir = DATA_DIR, work_dir = joinpath(dir, "mf6"), output_dir = joinpath(dir, "output"),
           start_date = start, end_date = stop, period = "month", initial_head_file = "",
           write_monthly_output = true, kwargs...)
end

"""
Indus-specific orientation fingerprints, on top of the generic ones in
`audit_static_inputs` (verified 2026-09-28 and 2026-10-16; examples/indus/README.md).
"""
function indus_fingerprints(s, p)
    sukkur = CartesianIndex(argmin(abs.(s.lon .- 68.87)), argmin(abs.(s.lat .- 27.7)))
    window = CartesianIndices((sukkur[1]-2:sukkur[1]+2, sukkur[2]-2:sukkur[2]+2))
    rowmean_area = vec(sum(s.cell_area; dims = 1)) ./ size(s.cell_area, 1)
    return [
        "grid is 180 rows x 204 columns" => size(s) == (204, 180),
        "cell area decreases northward" => issorted(rowmean_area; rev = true),
        "bankfull discharge > 5000 m3/s near Sukkur" => maximum(s.qbank[window]) > 5000,
        "Indus main stem in one catchment" => let ids = s.catchment[s.qbank .> 1000]
            maximum(count(==(u), ids) for u in unique(ids)) / length(ids) > 0.9
        end,
        "DRN on at least 95% of active land cells" =>
            length(p.drain.cells) >= 0.95 * count(s.landmask .& (s.idomain .> 0)),
    ]
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

@testset "ModflowCoupling: Indus example" begin
    tmp = mktempdir()
    cfg = test_config(tmp)
    s = read_static_inputs(cfg)
    p = build_parameters(s, cfg)

    @testset "Case 0: static audit" begin
        for (name, ok) in [audit_static_inputs(s, p, cfg); indus_fingerprints(s, p)]
            @test ok
        end
        @test 2 * length(p.chd_cells) == 1410          # reference_cases.md regression count
        @test length(p.river.cells) == 15165
        @test all(p.chd_head .== s.steady_head[p.chd_cells, :])
    end

    @testset "parameter formulas" begin
        @test MC.permeability_to_conductivity(-12.0) ≈ 1000 * 1e-12 * (9.81 / 0.001) * 86400
        r = p.river
        i = findfirst(r.is_river)
        @test MC.river_stage(r, i, -1.0) == r.bottom[i]            # no river infiltration
        @test MC.river_stage(r, i, 0.0) == r.bottom_bankfull[i]    # depth < 0.01 m -> 0
        q = 1000.0
        depth = (0.045q / (r.width[i] * sqrt(r.slope[i])))^(3 / 5)
        @test MC.river_stage(r, i, q) ≈ r.bottom_bankfull[i] + depth
        j = findfirst(!, r.is_river)
        @test MC.river_stage(r, j, 50.0) == r.head_ini[j]
        c = p.evt_cells[1]
        @test p.evt_rate[1] ≈ s.cprate[c] * s.cell_area[c] / (s.delr * s.delc)
    end

    @testset "Case 1/2: continuous run, BMI exchange" begin
        m = GroundwaterModel(cfg)
        @test BMI.get_current_time(m) == 0.0
        @test BMI.get_end_time(m) == 90.0
        @test BMI.get_time_step(m) == 31.0
        @test sort(BMI.get_input_var_names(m)) == ["groundwater_abstraction", "groundwater_recharge", "river_discharge"]
        riv, evt, rch, drn = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
        heads = Array{Float64,3}[]
        for frac in FORCING
            force!(m, frac)
            BMI.update(m)
            push!(riv, package_flows(m, "RIV")); push!(evt, package_flows(m, "EVT")); push!(rch, package_flows(m, "RCH"))
            push!(drn, package_flows(m, "DRN"))
            push!(heads, copy(m.head))
            # MODFLOW -> VIC conversions (coupling_contract.md)
            bf = reshape(BMI.get_value(m, "baseflow_to_routing", zeros(length(m.baseflow))), size(m.baseflow))
            @test all(isnan, bf[.!s.landmask])
            expected = zeros(size(bf))
            for (k, c) in enumerate(p.river.cells); expected[c] -= riv[end][k] / 86400; end
            for (k, c) in enumerate(p.drain.cells); expected[c] -= drn[end][k] / 86400; end
            @test all(bf[s.landmask] .≈ expected[s.landmask])
            @test all(<=(0), drn[end])                     # drains only remove water
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
        lst = joinpath(cfg.work_dir, lowercase(cfg.model_name) * ".lst")
        for (pkg, flows) in (("RIV", riv), ("EVT", evt), ("RCH", rch), ("DRN", drn))
            for (k, (qin, qout)) in enumerate(listing_rates(lst, pkg))
                @test sum(max.(flows[k], 0)) ≈ qin rtol = 1e-8 atol = 1e-3
                @test -sum(min.(flows[k], 0)) ≈ qout rtol = 1e-8 atol = 1e-3
            end
        end

        # Monthly output file
        NCDataset(joinpath(cfg.output_dir, "modflow_monthly.nc")) do ds
            @test length(ds["time"]) == 3
            @test isequal(coalesce.(ds["groundwater_head"][:, :, :, 3], NaN), heads[3])
        end

        @testset "restart after January equals continuous run" begin
            restart = test_config(joinpath(tmp, "restart"); start = Date(1979, 2, 1),
                                  initial_head_file = joinpath(cfg.output_dir, "head_19790131.nc"))
            m2 = GroundwaterModel(restart)
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
        include(joinpath(INDUS_DIR, "make_synthetic_forcing.jl"))
        forcing = make_synthetic_forcing(joinpath(tmp, "forcing", "vic.nc"), 1979:1979; config = cfg)
        cfgf = test_config(joinpath(tmp, "forcing"); stop = Date(1979, 2, 28), forcing_file = forcing)
        m = GroundwaterModel(cfgf)
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
        m2 = GroundwaterModel(cfgf)
        @test m.mf === nothing
        BMI.update(m2)
        BMI.finalize(m2)
    end

    @testset "input validation" begin
        cfgw = test_config(joinpath(tmp, "wel"); stop = Date(1979, 1, 31), pumping = true)
        m = GroundwaterModel(cfgw)
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

    @testset "config" begin
        # Relative paths resolve against the config's directory
        @test INDUS.data_dir == normpath(INDUS_DIR, "..", "..", "..", "data")
        @test INDUS.initial_head_file == joinpath(INDUS.data_dir, "steady_head_mgv_coupled.nc")
        @test INDUS.period == "day" && INDUS.model_name == "INDUS"
        @test INDUS.inputs["dem"].lat_order == :south_first
        @test INDUS.inputs["aqdepth"] == InputSpec("aqdepth_ini_indus.nc")
        # The synthetic config describes the same model, monthly from the global steady state
        syn = load_config(joinpath(INDUS_DIR, "synthetic.toml"))
        @test syn.inputs == INDUS.inputs && syn.data_dir == INDUS.data_dir
        @test syn.work_dir == joinpath(INDUS_DIR, "run", "mf6")
        @test syn.forcing_file == joinpath(INDUS_DIR, "vic_synthetic_monthly.nc")
        @test syn.period == "month" && isempty(syn.initial_head_file)

        dir = mktempdir()
        write_toml(name, text) = (path = joinpath(dir, name); write(path, text); path)
        base = """
            [paths]
            data_dir = "data"
            work_dir = "mf6"
            output_dir = "out"
            [time]
            start_date = 1979-01-01
            end_date = 1979-01-31
            [inputs]
            ibound = "ibound.nc"
            """
        @test_throws r"missing inputs: landmask" load_config(write_toml("base.toml", base))
        @test_throws r"unknown inputs landmsk" load_config(write_toml("typo.toml", base * "landmsk = \"x.nc\"\n"))
        @test_throws r"lat_order must be" load_config(write_toml("bad.toml", base * "dem = { file = \"x.nc\", lat_order = \"up\" }\n"))
        @test_throws r"unknown keys extends" load_config(write_toml("ext.toml", "extends = \"base.toml\"\n" * base))
        nohead = Config(cfg; inputs = filter(kv -> kv[1] != "steady_head", cfg.inputs))
        @test_throws r"initial_head_file or inputs.steady_head" MC.validate(nohead)
        nofoc = Config(cfg; foc = false, inputs = filter(kv -> !startswith(kv[1], "cp"), cfg.inputs))
        @test MC.validate(nofoc) === nofoc
    end

    @testset "sub-domain from larger maps" begin
        # A window of the Indus grid, with ibound stored north-first: every other
        # input is cut from the full-domain files by coordinate, including the
        # files with wrong latitude labels.
        dir = mktempdir()
        ii, jj = 61:120, 91:140
        NCDataset(joinpath(dir, "ibound_window.nc"), "c") do ds
            defVar(ds, "lon", s.lon[ii], ("lon",))
            defVar(ds, "lat", reverse(s.lat[jj]), ("lat",))
            defVar(ds, "ibound", reverse(s.idomain[ii, jj]; dims = 2), ("lon", "lat"))
        end
        inputs = merge(cfg.inputs, Dict("ibound" => InputSpec(joinpath(dir, "ibound_window.nc"))))
        w = read_static_inputs(Config(cfg; inputs))
        @test w.lon == s.lon[ii] && w.lat ≈ s.lat[jj]
        @test w.idomain == s.idomain[ii, jj]
        for f in (:landmask, :dem, :aqdepth, :cell_area, :qbank, :z0_floodplain, :cpdistinct, :catchment, :steady_head)
            full = getfield(s, f)
            @test isequal(getfield(w, f), ndims(full) == 3 ? full[ii, jj, :] : full[ii, jj])
        end
        @test w.delr ≈ s.delr && w.delc ≈ s.delc
    end
end
