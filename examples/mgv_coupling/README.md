# MODFLOW → mGV coupling example

`couple_mgv.jl` runs IndusModflow and [mGV](https://github.com/SWIFT-WUR/mGV)
in one Julia process and passes MODFLOW's monthly feedback to mGV through BMI.
It replaces the file exchange of the Python plugin for this direction.

| IndusModflow output | mGV BMI variable | Python plugin equivalent |
| --- | --- | --- |
| `baseflow_to_routing` (m3 s-1) | `channel_water__area~cross-sectional` (m2): `+= baseflow · dt / channel__length` before every daily mGV step | `discharge_mf` forcing for the VIC-WUR routing plugin |
| `capillary_rise` (mm month-1) | `soil__moisture_layer_N`, bottom layer (mm); excess moved up a layer at a time against `soil__moisture~maximum_layer_N` | `vic_runner.update_statefile` on `STATE_SOIL_MOISTURE` |

MODFLOW month M feeds mGV month M+1. Baseflow is held for the whole month.
Before every mGV time step, one step of baseflow volume (`baseflow · dt`) is
added to the channel storage. mGV stores a channel as its wetted
cross-section `area`, so storage is `area · length`. The kinematic wave then
releases the water downstream. Capillary rise is spread evenly over the days
of the month and added before every mGV time step. (The Python plugin added it
once, at the start of the month.)

The VIC → MODFLOW direction is not coupled: MODFLOW reads the offline
`forcing_file` from its config.

## Setup

mGV needs its `bmi-routing-vars` branch (tested at `7e4e4a6`), which exposes
all the variables above through BMI. They don't have to be listed in the
`[API]` section of the mGV config: that list only sets which names mGV
advertises.

```bash
cd /path/to/mGV
git switch bmi-routing-vars

cd /path/to/IndusModflow
julia --project=examples/mgv_coupling -e 'using Pkg; Pkg.develop([PackageSpec(path="."), PackageSpec(path="/path/to/mGV")])'
julia --project=examples/mgv_coupling examples/mgv_coupling/couple_mgv.jl \
    /path/to/mGV/configs/indus_config.toml examples/synthetic/synthetic.toml
```

Both models must start on the same date. The grids are matched by cell-centre
coordinates (within 0.001°). mGV cells outside the MODFLOW grid get zero
feedback.

## Caveats

- mGV still routes its own `subsurface_runoff`. In a two-way coupling, VIC
  baseflow and MODFLOW baseflow must not both reach the river.
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
  When channels run dry, only part of the loss is applied. On the Mekong
  case, −1 m3/s per cell for 30 days applied 1.2e9 of the requested
  2.4e9 m3. `add_baseflow!` returns the volume it applied.
- Cells without a channel length (NaN in the routing input) take no baseflow.
- Capillary rise that does not fit in the top layer is dropped, as in the
  Python plugin. `add_capillary_rise!` returns that volume.
- The script calls `BMI.update` once per day rather than `BMI.update_until`,
  because the baseflow and capillary rise are added before every step.
