#= Steady-state groundwater heads for the mGV <-> MODFLOW coupling.

How to run:
    julia --project=examples/mgv_coupling examples/mgv_coupling/spinup_steady_state.jl \
        <mgv_config.toml> <modflow_config.toml> <output_head.nc> [iterations] [skip_years]

The MODFLOW config must run daily periods (`period = "day"`); its
`initial_head_file` is the first starting point (the global steady state when
empty). The mGV config sets the years of each coupled run.

Each iteration:
  1. runs mGV and MODFLOW coupled, exchanging daily, over the mGV config's years;
  2. averages the recharge MODFLOW applied and the river discharge it received,
     leaving out the first `skip_years` (default 1) while mGV's soil settles;
  3. solves MODFLOW in steady state with those means, starting from the
     current heads, and writes the heads to `output_head.nc`.
The next iteration starts the coupled run from those heads. Point
`initial_head_file` of later runs at `output_head.nc` to start from them.

The steady state uses the mean river stage and lets rivers lose water without
checking mGV's channels, so its heads include the water the coupling creates.
=#
include(joinpath(@__DIR__, "coupling.jl"))
using Statistics, NCDatasets

mgv_config, mf_config, out_head = abspath.(ARGS[1:3])
iterations = parse(Int, get(ARGS, 4, "3"))
skip_years = parse(Int, get(ARGS, 5, "1"))

"`config` with some fields replaced."
with(config::Config; kw...) = Config(; (f => get(kw, f, getfield(config, f)) for f in fieldnames(Config))...)

# The steady-state solve starts far from its solution and can cycle where a river
# or EVT boundary switches on and off; damp the head changes more than the
# transient default (0.99) does.
const STEADY_SOLVER = IndusModflow.SolverSettings(outer_maximum = 2000, under_relaxation_gamma = 0.5)

"Mean forcing of the steady state, saved next to the head file."
function write_forcing(path, s, rch, q; attrib)
    NCDataset(path, "c"; attrib = attrib) do ds
        defVar(ds, "lon", s.lon, ("lon",)); defVar(ds, "lat", s.lat, ("lat",))
        defVar(ds, "groundwater_recharge", reshape(rch, size(s)), ("lon", "lat"); attrib = ["units" => "m d-1"])
        defVar(ds, "river_discharge", reshape(q, size(s)), ("lon", "lat"); attrib = ["units" => "m3 s-1"])
    end
end

base = load_config(mf_config)
base.period == "day" || error("$mf_config must run daily periods: set `period = \"day\"` under [time]")
start_head = base.initial_head_file
mean_abs(a, b) = mean(abs, filter(isfinite, a .- b))

for it in 1:iterations
    # 1. coupled run from the current heads
    gw = IndusGroundwater(with(base; initial_head_file = start_head, write_monthly_output = false))
    vic = BMI.initialize(mGV.Model, mgv_config)
    c = Coupling(gw, vic)
    head0 = copy(gw.head)
    first, last = current_date(gw), last_day(gw, vic)
    avg_from = first + Year(skip_years)
    avg_from < last || error("the mGV run is too short to skip $skip_years years")
    rch = zeros(size(c.vic_to_mf)); q = zeros(size(c.vic_to_mf)); n = 0
    t = @elapsed for day in first:Day(1):last
        discharge, recharge = step_day!(c, gw, vic)
        if day >= avg_from
            rch .+= max.(recharge, 0); q .+= discharge; n += 1
        end
    end
    drift = mean(filter(isfinite, gw.head[:, :, 1] .- head0[:, :, 1])) / ((last - first).value + 1) * 365.25
    years = ((last - first).value + 1) / 365.25
    @printf("iteration %d: coupled %s to %s in %.0f s; mean head drift %.3f m/year; water created %.1f km3/year\n",
            it, first, last, t, drift, created_water(c) / 1e9 / years)
    BMI.finalize(gw); BMI.finalize(vic)

    # 2-3. steady state with the mean forcing
    ss = IndusGroundwater(with(base; period = "steady", end_date = base.start_date,
                               initial_head_file = start_head, write_monthly_output = false,
                               work_dir = joinpath(base.work_dir, "steady"), solver = STEADY_SOLVER))
    BMI.set_value(ss, "groundwater_recharge", rch ./ n)
    BMI.set_value(ss, "river_discharge", q ./ n)
    prev = copy(ss.head)
    try
        BMI.update(ss)
    catch
        write_forcing(replace(out_head, r"\.nc$" => "") * "_forcing_failed.nc", ss.inputs, rch ./ n, q ./ n;
                      attrib = ["iteration" => it])
        rethrow()
    end
    change = mean_abs(ss.head[:, :, 1], prev[:, :, 1])
    attrib = ["source" => "spinup_steady_state.jl", "mgv_config" => mgv_config,
              "modflow_config" => mf_config, "iteration" => it,
              "mean_from" => string(avg_from), "mean_to" => string(last)]
    write_head(ss, out_head; long_name = "steady-state groundwater head, coupled mGV-MODFLOW iteration $it", attrib)
    write_forcing(replace(out_head, r"\.nc$" => "") * "_forcing.nc", ss.inputs, rch ./ n, q ./ n; attrib)
    @printf("iteration %d: steady state written to %s; mean |head change| %.2f m\n", it, out_head, change)
    BMI.finalize(ss)
    # the next iteration (and later runs) start from these heads
    global start_head = out_head
    flush(stdout)
end
