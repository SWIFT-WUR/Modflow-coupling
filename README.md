# ModflowCoupling.jl

MODFLOW 6 groundwater component for coupling with a land-surface model (mGV,
VIC) through BMI. It runs MODFLOW 6 in-process through
[ModflowInterface.jl](https://github.com/Deltares/ModflowInterface.jl) and
exposes a [BasicModelInterface.jl](https://github.com/Deltares/BasicModelInterface.jl)
API, so the land-surface model can exchange fields in memory. The domain and
its static inputs come from a TOML config: one basin, a sub-domain cut from
global maps, or the globe. It replaces the file-based Python VIC–MODFLOW plugin
(`vic_global_online_full.py`, `config_module.py`, `mf_run.py`), and the formulas
follow its coupling contract (`coupling_contract.md` of the Indus VIC–MODFLOW
handoff). Neither the plugin nor the handoff is part of this repository.

Examples:

- [`examples/indus/`](examples/indus/README.md): the Indus basin at 5 arcmin, with
  its config, synthetic forcing and the notes on its input data. The test suite
  runs this case.
- [`examples/mgv_coupling/`](examples/mgv_coupling/README.md): daily two-way
  coupling with [mGV](https://github.com/SWIFT-WUR/mGV), and a steady-state spin-up.

## Quick start

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'   # needs the Indus data, see examples/indus
julia --project=. examples/indus/make_synthetic_forcing.jl        # writes examples/indus/vic_synthetic_monthly.nc (not committed)
julia --project=. scripts/run_offline.jl examples/indus/synthetic.toml
```

```julia
using ModflowCoupling
import BasicModelInterface as BMI

gw = BMI.initialize(GroundwaterModel, "examples/indus/indus.toml")
BMI.set_value(gw, "groundwater_recharge", recharge)   # m/day on the MODFLOW grid
BMI.set_value(gw, "river_discharge", discharge)       # m3/s
BMI.update(gw)                                         # one stress period (a day in this config)
BMI.get_value_ptr(gw, "baseflow_to_routing")           # m3/s, for the next land-surface period
BMI.get_value_ptr(gw, "capillary_rise")                # mm per period, for the next land-surface period
BMI.finalize(gw)
```

`scripts/couple_vic_template.jl` shows the monthly loop with VIC-side
placeholders. `examples/mgv_coupling` is a working daily coupling.

## Configuration

A config is one TOML file that describes the whole MODFLOW model: domain,
inputs, time, packages and solver. Relative paths in `[paths]` are resolved
against the directory of the file, and input files against `data_dir`, so a
config moves together with its data. Unknown keys are an error.

```toml
model_name = "GWF"                 # MODFLOW 6 model name, at most 16 characters

[paths]
data_dir = "data"                  # where the [inputs] files are
work_dir = "run/mf6"               # MODFLOW 6 input and listing files
output_dir = "run/output"          # per-period output and restart heads
initial_head_file = ""             # optional restart head; default: inputs.steady_head
forcing_file = ""                  # optional offline forcing, {yyyy} and {mm} allowed

[time]
start_date = 1979-01-01
end_date = 1979-12-31
period = "month"                   # "month", "day" or "steady"

[coupling]
foc = true                         # EVT package + capillary-rise feedback
pumping = false                    # WEL package (human extension)

[output]
write_monthly = true
save_mf6_binary = false

[audit]
orientation = true                 # orientation fingerprints at every build

[forcing]                          # variable names in forcing_file
recharge_variable = "OUT_GWRECHARGE"   # mm/month
discharge_variable = "OUT_DISCHARGE"   # m3/s

[solver]                           # IMS settings; defaults follow the Python reference
outer_dvclose = 5.0

[inputs]
ibound = "ibound.nc"
dem = { file = "dem.nc", variable = "Band1", lat_order = "south_first" }
# ...
```

### Static inputs

Each entry of `[inputs]` is a file name, or a table with `file`, `variable`
and `lat_order`. Without `variable`, the file must hold exactly one data
variable. Variables have dimensions `(lat, lon)` in the file (`(lon, lat)` in
Julia order); coordinate names may be `lon`/`longitude`/`x` and
`lat`/`latitude`/`y`.

`ibound` defines the model grid. Every other input is selected at the ibound
cell centres by coordinate (within 0.001°), so it may cover a larger area at the
same resolution: for a small domain, point the inputs at global maps and give
only a smaller `ibound`. Latitude labels are trusted (`lat_order = "coordinates"`,
the default) unless an entry sets `lat_order = "south_first"` or
`"north_first"`: then the labels are known to be wrong, and the data are taken
to be stored in that order.

| Role | Required | Content |
| --- | --- | --- |
| `ibound` | yes | domain: 0 inactive, 1 active, 2 constant head (CHD); defines the grid |
| `landmask` | yes | land cells: finite and non-zero |
| `dem` | yes | land surface elevation (m), DIS top |
| `aqdepth` | yes | aquifer thickness (m) |
| `cell_area` | yes | true cell area (m²) |
| `ksat_log` | yes | log10 intrinsic permeability (m²) |
| `ksat_l1_conf_log`, `ksat_l2_conf_log` | no | log10 permeability of layer 1/2 where confined; `ksat_log` elsewhere |
| `spe_yi` | yes | specific yield |
| `qbank` | yes | bankfull discharge (m³/s) |
| `riv_slope` | yes | river slope |
| `min_dem` | yes | minimum elevation in the cell (m) |
| `z0_floodplain` | yes | floodplain elevation (m): drain level and river-bottom reference |
| `kq3` | yes | groundwater recession constant (1/day), DRN conductance |
| `catchment` | yes | catchment ids: a drain is not below the lowest `z0_floodplain` of its catchment |
| `cpsurface`, `cprate`, `cpdistinct` | with `foc` | EVT surface (m), maximum rate (m/day), extinction depth (m) |
| `steady_head` | no | two-layer head (m), `(lon, lat, layer)`: initial head and CHD head; needed unless `initial_head_file` is set and there are no CHD cells |

The DIS grid is uniform with a spacing of degrees × 111 km, as in the Python
reference. RCH and EVT rates are scaled by true cell area / DIS cell area, so
volumes follow `cell_area`.

`audit_static_inputs` runs on every model build and stops it when a check
fails. It checks the domain, layer geometry and package arrays. With
`orientation = true` under `[audit]` it also checks that the cell area follows
cos(latitude), that the DEM correlates with `steady_head` (r > 0.8), and that
bankfull discharge > 10 m³/s lies only on land. Add the fingerprints that are
specific to a dataset to its own tests, as `test/runtests.jl` does for the Indus.

### Offline forcing

Set `paths.forcing_file` to a NetCDF of monthly land-surface output (recharge
in mm/month, discharge in m3/s, dimensions lat/lon/time; VIC names by default,
see `[forcing]`). `{yyyy}` and `{mm}` in the path select one file per month.
Each `update` then reads the month itself; a field set through BMI since the
previous update takes precedence for that month.

### Visualise with the BMI Controller

The `server/` environment serves the model over
[Remote BMI](https://github.com/eWaterCycle/remotebmi):

```bash
julia --project=server -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=server server/bmi_server.jl                  # port 50051 (BMI_PORT to change)
```

Then start the [BMI Controller](https://github.com/BSchilperoort/bmi-controller)
with `pnpm dev`. Its dev server proxies to `localhost:50051`. Initialize with
the absolute config path, e.g. `/path/to/ModflowCoupling/examples/indus/synthetic.toml`.
Show grid 0 variables such as `groundwater_depth`, `baseflow_to_routing` or
`capillary_rise`; grid 1 (`groundwater_head`, 2 layers) has no colour map in
the controller. NaN values (off-domain cells) arrive as JSON `null`.
Re-initializing from the controller finalizes the running MF6 instance first,
because libmf6 allows only one instance per process.

## BMI variables

Grid 0 is the ibound grid (uniform rectilinear, values at cell centres, origin
at the south-west cell centre). Grid 1 adds 2 layers. Values are flattened
longitude-fastest with latitude ascending (south-first). Use `align_to_grid` to
map land-surface fields by coordinate. Time is reported in CF units `days since
<start_date>`, and one `update` is one stress period: a calendar month by
default, as in the Python reference, or one day with `period = "day"` under
`[time]` in the config (`start_date` may then be any day). With daily periods
`capillary_rise` is in mm d-1, and the per-period output goes to
`modflow_daily.nc`, with restart heads at month ends. `period = "steady"` makes
one steady-state stress period: one `update` solves for the heads in balance
with the inputs set through BMI, and `write_head(model, path)` saves them in the
form `initial_head_file` reads (`examples/mgv_coupling/spinup_steady_state.jl`
uses this).

| Variable | Direction | Units | Notes |
| --- | --- | --- | --- |
| `groundwater_recharge` | in | m d-1 | Mean over the period; non-finite values become 0 |
| `river_discharge` | in | m3 s-1 | Mean over the period; negative values set stage = river bottom |
| `groundwater_abstraction` | in | m3 d-1 | Needs `pumping = true`; ≥ 0 and only on `idomain == 1` |
| `groundwater_head` | out | m | End of period, both layers (grid 1); NaN on inactive cells |
| `groundwater_head_layer_1`, `_layer_2` | out | m | The same heads per layer, on grid 0 |
| `groundwater_depth` | out | m | Land surface (DIS top) minus layer-1 head; NaN off the landmask |
| `baseflow_to_routing` | out | m3 s-1 | `-(RIV + DRN flow)/86400`; NaN outside the landmask |
| `capillary_rise` | out | mm month-1 (mm d-1 with daily periods) | `-(EVT flow)·days/true area·1000`; FOC only |

Outputs of period M are meant for land-surface period M+1.

## Design

- **One continuous MF6 simulation.** Each calendar month (or day) is one stress period
  of one time step. The input files are written once, and each period
  `prepare_time_step`, then the stresses are set through BMI, then a
  `solve` loop runs until convergence (non-convergence is an error). Heads stay
  in memory in Float64. The Python reference rebuilt MF6 every month and
  passed heads through Float32 NetCDF.
- **BMI-written packages.** RCH (all active cells), RIV (static cell set; only
  the stage changes) and WEL hold all candidate records in PERIOD 1. MF6 6.3
  multiplies RCH rates by the cell area when it reads the PERIOD block.
  `calibrate_bmi!` detects this from a placeholder rate, so rates written
  through BMI are scaled the same way.
- **DRN.** A drain in every active land cell at the floodplain elevation `Z0`
  (not below the lowest `Z0` of the catchment), with conductance
  `KQ3 · storage coefficient · true cell area`, so groundwater above the
  floodplain drains like a linear reservoir. Its outflow is part of
  `baseflow_to_routing`. Records with a non-finite elevation or KQ3 are left out.
- **River bottom from `Z0`.** The river bottom is the floodplain elevation
  `Z0` minus the bankfull depth, not `min_dem` minus the bankfull depth as in
  the Python reference, so rivers do not lie above their own drains
  (`examples/indus/README.md`).
- **Restart.** Each month writes `output/head_YYYYMMDD.nc` (variable `gwl`). Set
  `paths.initial_head_file` to continue from it. The test suite checks that a
  restart after January reproduces the continuous February–March run
  (< 1e-6 m).
- **Per-period diagnostics** go to `output/modflow_monthly.nc` (`modflow_daily.nc`
  with daily periods).
- **MODFLOW 6 version.** 6.3.0, the version Modflow6_jll provides.

The Python reference behaviours below are kept on purpose: the IMS settings
(`outer_dvclose = 5 m`), NaN → 0 in restart heads, the steady head as the
first IC and CHD head, and scaling RCH/EVT by true area / DIS area.

## Not yet done

- Human case: capped pumping (`pumpingCapacity_...nc`, reference Case 3). WEL
  and its validation exist; the capacity capping does not.
- Longitude wrap-around: inputs are selected by coordinate, so a domain that
  crosses the edge of the input maps' longitude range (e.g. 180° for −180–180
  maps) cannot be cut from them.
- The VIC side: routing forcing and the soil-moisture update from
  `capillary_rise` (`vic_runner.update_statefile`) belong in the Julia VIC.
