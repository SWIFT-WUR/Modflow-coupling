"Views into MODFLOW 6 memory, valid between initialize and finalize."
struct Mf6Pointers
    head::Vector{Float64}                        # X, reduced nodes
    nodeuser::Vector{Int}                        # reduced node → user node
    rch_bound::Matrix{Float64}                   # (1, nrch): recharge rate
    riv_bound::Matrix{Float64}                   # (3, nriv): stage, cond, rbot
    riv_simvals::Vector{Float64}                 # m3/day, positive into the aquifer
    evt_simvals::Union{Vector{Float64},Nothing}  # m3/day, negative out of the aquifer
    wel_bound::Union{Matrix{Float64},Nothing}    # (1, nwel): q
    mxiter::Int
end

"""
    IndusGroundwater

MODFLOW 6 groundwater model of the Indus basin, advanced one stress period (a
calendar month, or a day with `period = "day"`) per `BMI.update`. Set the
period's VIC fields before calling `update`; afterwards the outputs describe the
end of that period and are meant for the *next* VIC period.
"""
mutable struct IndusGroundwater
    config::Config
    inputs::StaticInputs
    params::StaticParameters
    periods::Vector{Date}
    mf::Union{ModflowModel,Nothing}
    ptr::Union{Mf6Pointers,Nothing}
    # inputs for the coming stress period
    recharge::Matrix{Float64}       # m/day
    discharge::Matrix{Float64}      # m3/s
    abstraction::Matrix{Float64}    # m3/day
    # outputs of the last completed stress period
    head::Array{Float64,3}          # m
    baseflow::Matrix{Float64}       # m3/s
    capillary::Matrix{Float64}      # mm per stress period
    completed::Int
    rch_scale::Vector{Float64}      # MF6 BOUND value per unit RCH rate (1 or the DIS cell area)
    inputs_set::Set{Symbol}         # inputs set since the last update (override `forcing_file`)
end

"The model owning the process-wide libmf6 instance."
const ACTIVE_MODEL = Ref{Any}(nothing)

"""
    IndusGroundwater(config)

Read inputs, build parameters, write the MF6 simulation and initialize libmf6.
Only one MODFLOW 6 instance can be active per Julia process; a model that is
still active is finalized first.
"""
function IndusGroundwater(config::Config)
    ACTIVE_MODEL[] === nothing || finalize!(ACTIVE_MODEL[])
    inputs = read_static_inputs(config)
    params = build_parameters(inputs, config)
    failed = [name for (name, ok) in audit_static_inputs(inputs, params) if !ok]
    isempty(failed) || error("static input audit failed:\n  " * join(failed, "\n  "))
    periods = stress_periods(config)
    write_simulation(config, inputs, params, periods)

    grid = size(inputs)
    head = fill(NaN, grid..., 2)
    active = inputs.idomain .> 0
    for l in 1:2
        view(head, :, :, l)[active] .= view(params.initial_head, :, :, l)[active]
    end
    model = IndusGroundwater(
        config, inputs, params, periods, nothing, nothing,
        zeros(grid), zeros(grid), zeros(grid),
        head, fill(NaN, grid), zeros(grid), 0, Float64[], Set{Symbol}(),
    )
    initialize_mf6!(model)
    ACTIVE_MODEL[] = model
    return model
end

function mf6_pointer(mf, component, var; sub = "")
    return BMI.get_value_ptr(mf, MF.get_var_address(mf, var, component; subcomponent_name = sub))
end

function initialize_mf6!(m::IndusGroundwater)
    mf = ModflowModel(m.config.work_dir)
    BMI.initialize(mf)
    m.mf = mf
    name = m.config.model_name
    ptr(var, sub = "") = mf6_pointer(mf, name, var; sub)

    head = ptr("X")
    nuser = length(m.inputs.idomain) * 2
    nodeuser = Int.(vec(ptr("NODEUSER", "DIS")))
    if length(nodeuser) != length(head)       # no inactive cells: MF6 keeps a dummy array
        length(head) == nuser || error("unexpected MF6 node count $(length(head))")
        nodeuser = collect(1:nuser)
    end
    p = m.params
    evt_simvals = m.config.foc ? vec(ptr("SIMVALS", "EVT")) : nothing
    wel_bound = m.config.pumping ? ptr("BOUND", "WEL") : nothing
    m.ptr = Mf6Pointers(
        vec(head), nodeuser,
        ptr("BOUND", "RCH"), ptr("BOUND", "RIV"),
        vec(ptr("SIMVALS", "RIV")), evt_simvals, wel_bound,
        Int(only(mf6_pointer(mf, "SLN_1", "MXITER"))),
    )

    return m
end

"""
Run after the first `prepare_time_step`, when MF6 has read PERIOD 1: check that
the BMI boundary arrays follow the record order of `write_simulation`, and
record how MF6 stored the placeholder RCH rate of 1.0 (MF6 6.3 multiplies it by
the cell area on reading, so a rate set through BMI must be scaled the same way).
"""
function calibrate_bmi!(m::IndusGroundwater)
    p, ptr = m.params, m.ptr
    lin = LinearIndices(m.inputs.idomain)
    function check(pkg, cells)
        nodelist = vec(mf6_pointer(m.mf, m.config.model_name, "NODELIST"; sub = pkg))
        users = [ptr.nodeuser[n] for n in nodelist[1:length(cells)]]
        users == [lin[c] for c in cells] || error("$pkg record order differs from the written input")
    end
    check("RCH", p.rch_cells)
    check("RIV", p.river.cells)
    m.config.foc && check("EVT", p.evt_cells)
    m.config.pumping && check("WEL", p.wel_cells)
    m.rch_scale = ptr.rch_bound[1, :]
    all(x -> x == 1.0 || isapprox(x, m.inputs.delr * m.inputs.delc; rtol = 1e-12), m.rch_scale) ||
        error("unexpected RCH BOUND after reading period 1")
    return m
end

"Date of the first day of the next stress period to run."
current_date(m::IndusGroundwater) =
    m.completed < length(m.periods) ? m.periods[m.completed + 1] : m.periods[end] + period_step(m.config)

function set_recharge!(m::IndusGroundwater, recharge::AbstractArray)
    m.recharge .= reshape(recharge, size(m.recharge))
    replace!(x -> isfinite(x) ? x : 0.0, m.recharge)
    push!(m.inputs_set, :recharge)
    return m
end

function set_discharge!(m::IndusGroundwater, discharge::AbstractArray)
    m.discharge .= reshape(discharge, size(m.discharge))
    replace!(x -> isfinite(x) ? x : 0.0, m.discharge)
    push!(m.inputs_set, :discharge)
    return m
end

function set_abstraction!(m::IndusGroundwater, abstraction::AbstractArray)
    a = reshape(abstraction, size(m.abstraction))
    all(isfinite, a) || error("groundwater abstraction contains non-finite values")
    all(>=(0), a) || error("groundwater abstraction must be a non-negative extraction rate")
    any((a .> 0) .& (m.inputs.idomain .!= 1)) &&
        error("positive groundwater abstraction is only allowed on idomain == 1 cells")
    any(a .> 0) && !m.config.pumping && error("abstraction given but pumping is disabled")
    m.abstraction .= a
    push!(m.inputs_set, :abstraction)
    return m
end

"Write the period's RCH, RIV and WEL values into MF6 memory."
function apply_stresses!(m::IndusGroundwater)
    p, ptr = m.params, m.ptr
    for (i, c) in enumerate(p.rch_cells)
        ptr.rch_bound[1, i] = max(0.0, m.recharge[c] * p.area_factor[c]) * m.rch_scale[i]
    end
    r = p.river
    for (i, c) in enumerate(r.cells)
        ptr.riv_bound[1, i] = river_stage(r, i, m.discharge[c])
    end
    if ptr.wel_bound !== nothing
        for (i, c) in enumerate(p.wel_cells)
            ptr.wel_bound[1, i] = -m.abstraction[c]
        end
    end
end

"Read heads and exchange fluxes of the completed stress period."
function collect_outputs!(m::IndusGroundwater)
    p, ptr, s = m.params, m.ptr, m.inputs
    fill!(m.head, NaN)
    for (n, u) in zip(eachindex(ptr.head), ptr.nodeuser)
        m.head[u] = ptr.head[n]
    end
    for l in 1:2
        view(m.head, :, :, l)[s.idomain .<= 0] .= NaN
    end

    fill!(m.baseflow, 0.0)
    for (i, c) in enumerate(p.river.cells)
        m.baseflow[c] = -ptr.riv_simvals[i] / SECONDS_PER_DAY
    end
    replace!(x -> abs(x) > 1e8 ? 0.0 : x, m.baseflow)
    m.baseflow[.!s.landmask] .= NaN

    fill!(m.capillary, 0.0)
    if ptr.evt_simvals !== nothing
        days = period_days(m.config, m.periods[m.completed])
        for (i, c) in enumerate(p.evt_cells)
            m.capillary[c] = -ptr.evt_simvals[i] * days / s.cell_area[c] * 1000.0
        end
    end
    return m
end

"Advance one stress period."
function advance!(m::IndusGroundwater)
    m.completed < length(m.periods) || error("simulation already reached $(m.config.end_date)")
    date = m.periods[m.completed + 1]
    if !isempty(m.config.forcing_file)
        recharge, discharge = read_forcing(m.config, m.inputs, date)
        :recharge in m.inputs_set || set_recharge!(m, recharge)
        :discharge in m.inputs_set || set_discharge!(m, discharge)
    end
    empty!(m.inputs_set)
    mf = m.mf
    MF.prepare_time_step(mf, 0.0)
    m.completed == 0 && calibrate_bmi!(m)
    apply_stresses!(m)
    MF.prepare_solve(mf)
    converged = false
    for _ in 1:m.ptr.mxiter
        converged = MF.solve(mf)
        converged && break
    end
    MF.finalize_solve(mf)
    MF.finalize_time_step(mf)
    converged || error("MODFLOW 6 did not converge for the period starting $date")
    m.completed += 1
    collect_outputs!(m)
    m.config.write_monthly_output && write_monthly_output(m)
    return m
end

"Raw MF6 boundary flows (m3/day) of a package, in record order."
package_flows(m::IndusGroundwater, package::AbstractString) =
    copy(vec(mf6_pointer(m.mf, m.config.model_name, "SIMVALS"; sub = uppercase(package))))

function finalize!(m::IndusGroundwater)
    if m.mf !== nothing
        m.ptr = nothing
        BMI.finalize(m.mf)
        m.mf = nothing
    end
    ACTIVE_MODEL[] === m && (ACTIVE_MODEL[] = nothing)
    return m
end
