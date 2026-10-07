# MODFLOW → mGV coupling example

`couple_mgv.jl` runs IndusModflow and [mGV](https://github.com/SWIFT-WUR/mGV)
in one Julia process and passes MODFLOW's monthly feedback to mGV through BMI.
It replaces the file exchange of the Python plugin for this direction.

| IndusModflow output | mGV BMI variable | Python plugin equivalent |
| --- | --- | --- |
| `baseflow_to_routing` (m3 s-1) | `river_channel_area` (m2): `+= baseflow · dt / river_cell_length` before every daily mGV step | `discharge_mf` forcing for the VIC-WUR routing plugin |
| `capillary_rise` (mm month-1) | bottom layer of `moisture` (mm); excess moved up a layer at a time against `maximum_moisture` | `vic_runner.update_statefile` on `STATE_SOIL_MOISTURE` |

MODFLOW month M feeds mGV month M+1. Baseflow is held for the whole month.
Before every mGV time step, one step of baseflow volume (`baseflow · dt`) is
added to the channel storage. mGV stores a channel as its wetted
cross-section `area`, so storage is `area · length`. The kinematic wave then
releases the water downstream. Capillary rise is added once, at the start of
the month. `moisture` and `maximum_moisture` are read and written as whole
(lon, lat, layer) arrays, so the per-layer `_layer_N` names are not needed.

The VIC → MODFLOW direction is not coupled: MODFLOW reads the offline
`forcing_file` from its config.

## Setup

mGV needs its BMI branch (tested at `c005d91`) plus
`mgv_bmi_coupling_inputs.patch`. The patch only adds lookup-table entries, so
`river_channel_area`, `river_cell_length`, `river_discharge` and
`maximum_moisture` are reachable through BMI. It does not change mGV's physics:

```bash
cd /path/to/mGV
git switch 11-requested-feature-implement-bmi
git apply /path/to/IndusModflow/examples/mgv_coupling/mgv_bmi_coupling_inputs.patch

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
- The script adds the baseflow to channel storage, not to `river_discharge`.
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
  2.4e9 m3. The script prints the applied and requested volumes every month.
- Cells without a channel length (NaN in the routing input) take no baseflow.
- Capillary rise that does not fit in the top layer is dropped, as in the
  Python plugin. The script prints that volume each month.
- mGV's `BMI.update_until` does not work yet on the branch, so the script
  calls `BMI.update` once per day.
