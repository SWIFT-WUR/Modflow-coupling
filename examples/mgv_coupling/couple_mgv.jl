#= mGV <-> MODFLOW coupling example.

How to run:
    julia --project=examples/mgv_coupling examples/mgv_coupling/couple_mgv.jl \
        <mgv_config.toml> <modflow_config.toml>

The two models exchange every day. MODFLOW must run daily stress periods:
`period = "day"` in the [time] section of its config.

mGV -> MODFLOW (VIC -> MODFLOW in coupling_contract.md), the values of mGV day
D, used by MODFLOW day D:

  mGV `channel_water__volume_flow_rate` [m3/s] -> `river_discharge` [m3/s]
      Python: OUT_DISCHARGE. Sets the river stage of the RIV package.
  mGV `land_subsurface__runoff` [mm/day]      -> `groundwater_recharge` [m/day]
      Python: OUT_GWRECHARGE. The drainage out of the bottom soil layer. mGV
      must run with `external_groundwater_routing = true` in its [API] config,
      so this water reaches the river only through MODFLOW and not also through
      mGV's own routing.
  `groundwater_abstraction` is not coupled: mGV has no water-demand module, and
      the natural case of the contract has no pumping.

MODFLOW -> mGV, the two feedbacks of the Python plugin:

  baseflow_to_routing [m3/s]     -> mGV `channel_water__area~cross-sectional`
      Python: `discharge_mf`, read by VIC-WUR's routing plugin as MONTH forcing.
      Here: one day of baseflow volume is added to the channel storage before
      the next mGV step, and mGV's kinematic wave releases it downstream.
  capillary_rise      [mm/day]   -> mGV `soil__moisture_layer_N`, bottom layer (FOC only)
      Python: vic_runner.update_statefile, added to STATE_SOIL_MOISTURE layer 3
      once at the start of the month. Here: added before the next mGV step.
      Excess above `soil__moisture~maximum` moves up one layer at a time.

Timing: mGV runs day D with the feedback from MODFLOW day D-1, then MODFLOW runs
day D with mGV's day-D values, then its feedback is set for mGV day D+1.

A river loss (negative baseflow) can only be taken from water in mGV's channel.
The volume MODFLOW calculated but mGV could not supply enters the aquifer
without leaving the river; the script reports it at the end.
=#
include(joinpath(@__DIR__, "coupling.jl"))

mgv_config = abspath(get(ARGS, 1, ""))
mf_config = abspath(get(ARGS, 2, ""))
(isfile(mgv_config) && isfile(mf_config)) || error("usage: couple_mgv.jl <mgv_config.toml> <modflow_config.toml>")

gw = BMI.initialize(IndusGroundwater, mf_config)
vic = BMI.initialize(mGV.Model, mgv_config)
coupling = Coupling(gw, vic)

for day in current_date(gw):Day(1):last_day(gw, vic)
    step_day!(coupling, gw, vic)
end

@printf("baseflow requested %.4g m3, applied %.4g m3, not supplied by the channels %.4g m3\n",
        coupling.requested, coupling.applied, created_water(coupling))

BMI.finalize(gw)
BMI.finalize(vic)
