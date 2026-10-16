#=
Writes a MODFLOW 6 simulation covering all monthly stress periods. Package
settings follow mf_run.mfrun.run_modflow, plus the DRN package of the global
model (`get_drn_param`), which the Python Indus case left out. DRN is static.
Boundary packages whose values change each month (RCH, RIV, WEL) contain records for every candidate cell in period 1;
the values are overwritten through BMI before each time step, so later periods
need no PERIOD block. MF6 6.3 multiplies RCH rates by the cell area while
reading the PERIOD block, so values written through BMI must include that
factor (see `calibrate_bmi!`). The record order written here is the order of the BMI
`BOUND`/`SIMVALS` arrays.
=#

"Shortest round-trip decimal representation, readable by Fortran."
function fmt(x::Real)
    isfinite(x) || error("non-finite value $x in MODFLOW input")
    return x isa Integer ? string(x) : repr(Float64(x))
end

"Write a (ncol, nrow[, nlay]) array. MF6 reads each model row as a separate record."
function write_array(io, name, a::AbstractArray; per_line = 10)
    layered = ndims(a) == 3
    println(io, "  ", name, layered ? " LAYERED" : "")
    for layer in eachslice(reshape(a, size(a, 1), size(a, 2), :); dims = 3)
        println(io, "    INTERNAL FACTOR ", eltype(a) <: Integer ? "1" : "1.0")
        for row in eachcol(layer), chunk in Iterators.partition(row, per_line)
            println(io, "      ", join(fmt.(chunk), " "))
        end
    end
end

function write_block(f, io, name)
    println(io, "BEGIN ", name)
    f()
    println(io, "END ", name, "\n")
end

cellid(c::CartesianIndex{2}, layer = 1) = "$layer $(c[2]) $(c[1])"   # lay row col

function write_list_package(path, cells, values; options = String[], dims = String[], layers = nothing)
    open(path, "w") do io
        write_block(io, "OPTIONS") do
            foreach(o -> println(io, "  ", o), options)
        end
        write_block(io, "DIMENSIONS") do
            println(io, "  MAXBOUND ", length(cells))
            foreach(d -> println(io, "  ", d), dims)
        end
        write_block(io, "PERIOD 1") do
            for (i, c) in enumerate(cells)
                layer = layers === nothing ? 1 : layers[i]
                println(io, "  ", cellid(c, layer), " ", join(fmt.(values[i, :]), " "))
            end
        end
    end
end

"""
    write_simulation(config, inputs, params, periods)

Write the MODFLOW 6 input files into `config.work_dir`.
"""
function write_simulation(config::Config, s::StaticInputs, p::StaticParameters, periods::Vector{Date})
    dir = config.work_dir
    mkpath(dir)
    name = lowercase(config.model_name)
    ncol, nrow = size(s)
    sv = config.solver
    file(ext) = joinpath(dir, "$name.$ext")

    open(joinpath(dir, "mfsim.nam"), "w") do io
        write_block(() -> nothing, io, "OPTIONS")
        write_block(() -> println(io, "  TDIS6 $name.tdis"), io, "TIMING")
        write_block(() -> println(io, "  GWF6 $name.nam $(config.model_name)"), io, "MODELS")
        write_block(() -> nothing, io, "EXCHANGES")
        write_block(() -> println(io, "  IMS6 $name.ims $(config.model_name)"), io, "SOLUTIONGROUP 1")
    end

    open(file("tdis"), "w") do io
        write_block(() -> println(io, "  TIME_UNITS DAYS"), io, "OPTIONS")
        write_block(() -> println(io, "  NPER ", length(periods)), io, "DIMENSIONS")
        write_block(io, "PERIODDATA") do
            for d in periods
                println(io, "  ", fmt(Float64(period_days(config, d))), " 1 1.0")
            end
        end
    end

    open(file("ims"), "w") do io
        write_block(io, "OPTIONS") do
            println(io, "  PRINT_OPTION SUMMARY")
            println(io, "  COMPLEXITY COMPLEX")
        end
        write_block(io, "NONLINEAR") do
            println(io, "  OUTER_DVCLOSE ", fmt(sv.outer_dvclose))
            println(io, "  OUTER_MAXIMUM ", sv.outer_maximum)
            println(io, "  UNDER_RELAXATION SIMPLE")
            println(io, "  UNDER_RELAXATION_GAMMA ", fmt(sv.under_relaxation_gamma))
        end
        write_block(io, "LINEAR") do
            println(io, "  INNER_MAXIMUM ", sv.inner_maximum)
            println(io, "  INNER_DVCLOSE ", fmt(sv.inner_dvclose))
            println(io, "  INNER_RCLOSE ", fmt(sv.inner_rclose), " STRICT")
            println(io, "  LINEAR_ACCELERATION BICGSTAB")
            println(io, "  RELAXATION_FACTOR ", fmt(sv.relaxation_factor))
        end
    end

    packages = ["DIS6 $name.dis DIS", "NPF6 $name.npf NPF", "STO6 $name.sto STO", "IC6 $name.ic IC",
                "CHD6 $name.chd CHD", "RCH6 $name.rch RCH", "RIV6 $name.riv RIV", "DRN6 $name.drn DRN"]
    config.foc && push!(packages, "EVT6 $name.evt EVT")
    config.pumping && push!(packages, "WEL6 $name.wel WEL")
    push!(packages, "OC6 $name.oc OC")
    open(file("nam"), "w") do io
        write_block(() -> nothing, io, "OPTIONS")
        write_block(() -> foreach(l -> println(io, "  ", l), packages), io, "PACKAGES")
    end

    open(file("dis"), "w") do io
        write_block(io, "OPTIONS") do
            println(io, "  LENGTH_UNITS METERS")
            println(io, "  NOGRB")
        end
        write_block(io, "DIMENSIONS") do
            println(io, "  NLAY 2\n  NROW $nrow\n  NCOL $ncol")
        end
        write_block(io, "GRIDDATA") do
            println(io, "  DELR\n    CONSTANT ", fmt(s.delr))
            println(io, "  DELC\n    CONSTANT ", fmt(s.delc))
            write_array(io, "TOP", p.top)
            write_array(io, "BOTM", p.botm)
            write_array(io, "IDOMAIN", cat(s.idomain, s.idomain; dims = 3))
        end
    end

    open(file("npf"), "w") do io
        write_block(() -> nothing, io, "OPTIONS")
        write_block(io, "GRIDDATA") do
            println(io, "  ICELLTYPE\n    CONSTANT 0")
            write_array(io, "K", p.k)
            write_array(io, "K33", p.k33)
        end
    end

    open(file("sto"), "w") do io
        write_block(() -> println(io, "  STORAGECOEFFICIENT"), io, "OPTIONS")
        write_block(io, "GRIDDATA") do
            println(io, "  ICONVERT\n    CONSTANT 0")
            write_array(io, "SS", p.ss)
        end
        write_block(() -> println(io, config.period == "steady" ? "  STEADY-STATE" : "  TRANSIENT"), io, "PERIOD 1")
    end

    open(file("ic"), "w") do io
        write_block(() -> nothing, io, "OPTIONS")
        write_block(() -> write_array(io, "STRT", p.initial_head), io, "GRIDDATA")
    end

    save = config.save_mf6_binary ? ["SAVE_FLOWS"] : String[]
    nchd = length(p.chd_cells)
    write_list_package(file("chd"), vcat(p.chd_cells, p.chd_cells), reshape(vec(p.chd_head), :, 1);
                       options = save, layers = vcat(fill(1, nchd), fill(2, nchd)))
    # Placeholder rate 1.0: `calibrate_bmi!` reads back how MF6 stored it.
    write_list_package(file("rch"), p.rch_cells, ones(length(p.rch_cells), 1); options = save)
    r = p.river
    stage0 = [river_stage(r, i, 0.0) for i in eachindex(r.cells)]
    write_list_package(file("riv"), r.cells, hcat(stage0, r.conductance, r.bottom); options = save)
    d = p.drain
    write_list_package(file("drn"), d.cells, hcat(d.elevation, d.conductance); options = save)
    if config.foc
        write_list_package(file("evt"), p.evt_cells, hcat(p.evt_surface, p.evt_rate, p.evt_depth);
                           options = ["FIXED_CELL"; save], dims = ["NSEG 1"])
    end
    if config.pumping
        write_list_package(file("wel"), p.wel_cells, zeros(length(p.wel_cells), 1); options = save)
    end

    open(file("oc"), "w") do io
        write_block(io, "OPTIONS") do
            if config.save_mf6_binary
                println(io, "  HEAD FILEOUT $name.hds")
                println(io, "  BUDGET FILEOUT $name.cbc")
            end
        end
        write_block(io, "PERIOD 1") do
            println(io, "  PRINT BUDGET ALL")
            if config.save_mf6_binary
                println(io, "  SAVE HEAD ALL")
                println(io, "  SAVE BUDGET ALL")
            end
        end
    end
    return dir
end
