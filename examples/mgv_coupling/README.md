# mGV ↔ MODFLOW coupling example

Files:

- `couple_mgv.jl` runs the coupled models.
- `spinup_steady_state.jl` computes steady-state groundwater heads to start from
  (see "Starting heads" below).
- `coupling.jl` holds the functions both scripts share.
- `../indus/indus.toml` is the daily MODFLOW config for the Indus that starts
  from the stored steady-state heads. For another domain, write a config with
  that domain's `[inputs]` and `period = "day"` (see the main README).

`couple_mgv.jl` runs the ModflowCoupling MODFLOW model and [mGV](https://github.com/SWIFT-WUR/mGV)
in one Julia process and exchanges the fields of `coupling_contract.md` through
BMI, in both directions, every day. It replaces the monthly file exchange of
the Python plugin. MODFLOW runs daily stress periods for this: set
`period = "day"` in the `[time]` section of its config (the default, `"month"`,
matches the Python reference).

## mGV → MODFLOW

| mGV BMI variable | MODFLOW input | Python plugin equivalent |
| --- | --- | --- |
| `channel_water__volume_flow_rate` (m3 s-1) | `river_discharge` (m3 s-1) | `OUT_DISCHARGE` |
| `land_subsurface__runoff` (mm day-1) ÷ 1000 | `groundwater_recharge` (m day-1) | `OUT_GWRECHARGE` |

mGV's values of day D drive MODFLOW day D. `land_subsurface__runoff` is the
drainage out of the bottom soil layer. mGV must run with
`external_groundwater_routing = true` in the `[API]` section of its config, so
this water reaches the river only through MODFLOW (as baseflow) and not also
through mGV's own routing. The script stops with an error if it isn't set.
`groundwater_abstraction` is not coupled: mGV has no water-demand module, and
the natural case of the contract has no pumping.

## MODFLOW → mGV

| MODFLOW output | mGV BMI variable | Python plugin equivalent |
| --- | --- | --- |
| `baseflow_to_routing` (m3 s-1) | `channel_water__area~cross-sectional` (m2): `+= baseflow · dt / channel__length` before the next mGV step | `discharge_mf` forcing for the VIC-WUR routing plugin |
| `capillary_rise` (mm day-1) | `soil__moisture_layer_N`, bottom layer (mm); excess moved up a layer at a time against `soil__moisture~maximum_layer_N` | `vic_runner.update_statefile` on `STATE_SOIL_MOISTURE` |

MODFLOW day D feeds mGV day D+1. Before that mGV step, one day of baseflow
volume (`baseflow · dt`) is added to the channel storage. mGV stores a channel
as its wetted cross-section `area`, so storage is `area · length`. The
kinematic wave then releases the water downstream. Capillary rise is added to
the soil before the same step. (The Python plugin exchanged monthly, a month
later, and added capillary rise once at the start of the month.)

## Setup

mGV needs its `mgv-to-modflow-coupling` branch (tested at `34bccf3`), which
exposes all the variables above through BMI and has the
`external_groundwater_routing` setting. The variables don't have to be listed in
`variables` in the `[API]` section of the mGV config: that list only sets
which names mGV advertises.

```bash
cd /path/to/mGV
git switch mgv-to-modflow-coupling
# and set `external_groundwater_routing = true` under [API] in the mGV config

cd /path/to/ModflowCoupling
julia --project=examples/mgv_coupling -e 'using Pkg; Pkg.develop([PackageSpec(path="."), PackageSpec(path="/path/to/mGV")])'
julia --project=examples/mgv_coupling examples/mgv_coupling/couple_mgv.jl \
    /path/to/mGV/configs/indus_config.toml examples/indus/indus.toml
```

The MODFLOW config must have `period = "day"` under `[time]`; the script stops
with an error if it doesn't. The mGV config sets the years that run.

## Starting heads

The global steady-state heads that MODFLOW starts from by default are close to
balance with the coupled models, but not in it: from them the heads rise about
0.2 m/year at first. `spinup_steady_state.jl` finds the balance directly. Each
iteration runs the coupled models daily over the mGV config's years, solves
MODFLOW in steady state with the mean recharge and river discharge of that run
(leaving out the first year), and starts the next iteration from those heads:

```bash
julia --project=examples/mgv_coupling examples/mgv_coupling/spinup_steady_state.jl \
    /path/to/mGV/configs/indus_1975_1984.toml examples/indus/indus.toml \
    ../data/steady_head_mgv_coupled.nc 5
```

It writes the heads (variable `gwl`, the form `initial_head_file` reads) and
the mean forcing (`..._forcing.nc`). `indus.toml` points
`initial_head_file` at `../data/steady_head_mgv_coupled.nc`, which holds the
result for 1975–1984 forcing (GFDL-ESM4 historical) after five iterations,
starting from the global steady state, with the DRN package and river bottoms
at `Z0`:

| Iteration | Mean head change of the steady solve | Head drift in the next coupled run | Water created by the coupling |
| --- | --- | --- | --- |
| start | – | +0.19 m/year | 11.0 km3/year |
| 1 | 68 m | −0.004 m/year | 6.7 km3/year |
| 2 | 0.03 m | −0.005 m/year | 6.7 km3/year |
| 4 | < 0.01 m | −0.005 m/year | 6.7 km3/year |

Compared with the global steady state (14,476 variable-head land cells, top
layer), the heads correlate at r = 0.991. Half of the cells lie within 0.3 m of
the global head and 73% within 25 m. Below 500 m elevation (7,240 cells, the
plains) the median difference is 0.0 m (10th–90th percentile −7.6 to +5.4 m)
and the median depth to the water table is 20 m in both. The differences are in
the mountains: above 4,000 m the new heads are a median 38 m higher, up to
several hundred metres in places. The median depth to the water table is 121
m (global 122 m). The earlier spin-up without DRN, kept as
`../data/steady_head_mgv_coupled_noDRN.nc`, had a median depth of 7 m.

The steady-state solve uses stronger under-relaxation (0.5) than transient
periods, because from far away it cycles where a river or EVT boundary
switches on and off. Its rivers lose water at the mean stage without checking
mGV's channels, so the heads include the water the coupling creates.

Both models must start on the same date. The grids are matched by cell-centre
coordinates (within 0.001°), in both directions. mGV cells outside the MODFLOW
grid get zero feedback, and MODFLOW cells outside the mGV grid get zero recharge
and discharge (the script warns about them). Values set through BMI take
precedence over a `forcing_file` in the MODFLOW config.

## Caveats

- MODFLOW does not route the river. It uses mGV's discharge only to set the
  river stage of its RIV package, through its own channel geometry (width from
  bankfull discharge, the `riv_slope` input). That geometry differs from mGV's
  routing (width from upstream area, its own slopes), so MODFLOW's stage is not
  the depth of the water in mGV's channel.
- With `external_groundwater_routing = true`, mGV's rivers carry only surface
  runoff plus MODFLOW's baseflow from the day before. On the Mekong case in
  January, with monthly exchange, the largest monthly mean discharge was
  0.1 m3/s before that baseflow came back.
- The script adds the baseflow to channel storage, not to
  `channel_water__volume_flow_rate` (discharge).
  In mGV's kinematic wave, the stored discharge is only the previous outflow
  (`Q_old`). Raising it drains the channel, so discharge goes down: on the
  Mekong case, adding 10 m3/s per cell lowered the outlet from 140.5 to
  120.9 m3/s. Adding the same amount to storage raises discharge. With
  1 m3/s per channel cell for 30 days, the outlet went from 140.5 to
  458.6 m3/s. 339 channel cells drain to the outlet, so the fully settled
  rise would be 339 m3/s; after 30 days it was 318 m3/s.
- Negative baseflow (a losing river) cannot take the channel area below zero.
  MODFLOW calculates a day's river loss from its river stage without knowing
  how much water mGV's channel holds. Where the loss is more than the channel
  contains, the rest enters the aquifer without leaving the river: the
  coupling creates water. The script prints that volume at the end. Daily
  exchange roughly halves it: on the Indus case over 1979–1980, 126 km3 with
  daily exchange against 240 km3 with monthly exchange, both from the global
  steady-state heads and before the DRN package was added. Starting from the coupled steady state lowers it further,
  to about 7 km3/year, but does not remove it.
- Cells without a channel length (NaN in the routing input) take no baseflow.
- Capillary rise that does not fit in the top layer is dropped, as in the
  Python plugin. `add_capillary_rise!` returns that volume.
- The script calls `BMI.update` once per day rather than `BMI.update_until`,
  because the baseflow and capillary rise are added before every step.
