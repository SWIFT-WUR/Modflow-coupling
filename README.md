# IndusModflow.jl

MODFLOW 6 groundwater component for the Indus VIC–MODFLOW coupling. It runs
MODFLOW 6 in-process through [ModflowInterface.jl](https://github.com/Deltares/ModflowInterface.jl)
and exposes a [BasicModelInterface.jl](https://github.com/Deltares/BasicModelInterface.jl)
API, so a BMI-enabled Julia VIC can exchange fields in memory. It replaces the
file-based Python VIC–MODFLOW plugin (`vic_global_online_full.py`,
`config_module.py`, `mf_run.py`). The formulas follow the coupling contract of
the Indus VIC–MODFLOW reference handoff (`coupling_contract.md`). Neither the
plugin nor the handoff is part of this repository.

## Input data

The static inputs are not in this repository. The example configs and the tests
expect them in `../data` (a `data/` directory next to this checkout). The tests
also accept `INDUS_MODFLOW_DATA=/path/to/data`. The file names can be changed in
the `[files]` table of a config (see `DEFAULT_FILES` in `src/config.jl`).

| Role | File |
| --- | --- |
| idomain | `ibound_indus_normalized.nc` |
| grid spacing | `Clone_05min_indus.nc` |
| DIS top, aquifer depth | `dem_ini_indus.nc`, `aqdepth_ini_indus.nc` |
| true cell area | `cell_area_indus.nc` |
| permeability | `ksat_log_indus.nc`, `ksat_l1_conf_log_indus.nc`, `ksat_l2_conf_log_indus.nc` |
| storage | `spe_yi_inp_indus.nc` |
| land mask | `landmask_indus.nc` |
| rivers | `bankfull_discharge_indus.nc`, `riv_slope_indus.nc`, `min_dem_indus.nc`, `Z0_floodplain_indus_updated.nc` |
| capillary rise (EVT) | `capillary_rise_indus.nc` |
| initial and CHD head | `mf_ss_gwl_5minDRN.nc` (global steady state) |

## Quick start

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. scripts/make_synthetic_forcing.jl          # writes examples/synthetic/vic_synthetic_monthly.nc (not committed)
julia --project=. scripts/run_offline.jl examples/synthetic/synthetic.toml
```

```julia
using IndusModflow
import BasicModelInterface as BMI

gw = BMI.initialize(IndusGroundwater, "scripts/indus_natural_foc.toml")
BMI.set_value(gw, "groundwater_recharge", recharge)   # m/day on the MODFLOW grid
BMI.set_value(gw, "river_discharge", discharge)       # m3/s
BMI.update(gw)                                         # one calendar month
BMI.get_value_ptr(gw, "baseflow_to_routing")           # m3/s, for next VIC month
BMI.get_value_ptr(gw, "capillary_rise")                # mm/month, for next VIC month
BMI.finalize(gw)
```

`scripts/couple_vic_template.jl` shows the monthly loop. Its VIC-side calls are
placeholders.

### Offline forcing

Set `paths.forcing_file` to a NetCDF of monthly VIC output (`OUT_GWRECHARGE` in
mm/month, `OUT_DISCHARGE` in m3/s, dimensions lat/lon/time). `{yyyy}` and `{mm}`
in the path select one file per month. Each `update` then reads the month
itself; a field set through BMI since the previous update takes precedence for
that month. `scripts/make_synthetic_forcing.jl` writes illustrative forcing
(1979–1981; not a VIC simulation): a monsoon recharge peak in Jul–Sep that is
strongest in the north-east, and river discharge as a seasonal fraction of
bankfull discharge. The years are normal, wet and dry.

### Visualise with the BMI Controller

The `server/` environment serves the model over
[Remote BMI](https://github.com/eWaterCycle/remotebmi):

```bash
julia --project=server -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=server server/bmi_server.jl                  # port 50051 (BMI_PORT to change)
```

Then start the [BMI Controller](https://github.com/BSchilperoort/bmi-controller)
with `pnpm dev`. Its dev server proxies to `localhost:50051`. Initialize with
the absolute config path, e.g.
`/home/bart/data/modflow/IndusModflow/examples/synthetic/synthetic.toml`. Show
grid 0 variables such as `groundwater_depth`, `baseflow_to_routing` or
`capillary_rise`; grid 1 (`groundwater_head`, 2 layers) has no colour map in
the controller. NaN values (off-basin cells) arrive as JSON `null`.
Re-initializing from the controller finalizes the running MF6 instance first,
because libmf6 allows only one instance per process.

## BMI variables

Grid 0 is 180 lat × 204 lon (uniform rectilinear, values at cell centres,
origin at the south-west cell centre). Grid 1 adds 2 layers. Values are
flattened longitude-fastest with latitude ascending (south-first), matching
`ibound_indus_normalized.nc`. Use `align_to_grid` to map VIC fields by
coordinate. Time is reported in CF units `days since <start_date>`, and one
`update` is one month.

| Variable | Direction | Units | Notes |
| --- | --- | --- | --- |
| `groundwater_recharge` | in | m d-1 | Monthly mean; non-finite values become 0 |
| `river_discharge` | in | m3 s-1 | Monthly mean; negative values set stage = river bottom |
| `groundwater_abstraction` | in | m3 d-1 | Needs `pumping = true`; ≥ 0 and only on `idomain == 1` |
| `groundwater_head` | out | m | End of month, both layers (grid 1); NaN on inactive cells |
| `groundwater_head_layer_1`, `_layer_2` | out | m | The same heads per layer, on grid 0 |
| `groundwater_depth` | out | m | Land surface (DIS top) minus layer-1 head; NaN off the landmask |
| `baseflow_to_routing` | out | m3 s-1 | `-(RIV flow)/86400`; NaN outside the landmask |
| `capillary_rise` | out | mm month-1 | `-(EVT flow)·days/true area·1000`; FOC only |

Outputs of month M are meant for VIC month M+1.

## Design

- **One continuous MF6 simulation.** Each calendar month is one stress period
  of one time step. The input files are written once, and each month
  `prepare_time_step`, then the stresses are set through BMI, then a
  `solve` loop runs until convergence (non-convergence is an error). Heads stay
  in memory in Float64. The Python reference rebuilt MF6 every month and
  passed heads through Float32 NetCDF.
- **BMI-written packages.** RCH (all active cells), RIV (static cell set; only
  the stage changes) and WEL hold all candidate records in PERIOD 1. MF6 6.3
  multiplies RCH rates by the cell area when it reads the PERIOD block.
  `calibrate_bmi!` detects this from a placeholder rate, so rates written
  through BMI are scaled the same way.
- **Restart.** Each month writes `output/head_YYYYMMDD.nc` (variable `gwl`). Set
  `paths.initial_head_file` to continue from it. The test suite checks that a
  restart after January reproduces the continuous February–March run
  (< 1e-6 m).
- **Monthly diagnostics** go to `output/indus_modflow_monthly.nc`.

## Differences from the Python reference

1. **Grid orientation (important).** The static NetCDFs do not share one
   latitude storage order, and `dem_ini_indus.nc` and `landmask_indus.nc` have
   wrong latitude labels. The Python code combined raw `.values`. As a result,
   aquifer depth, cell area, K, storage, bankfull discharge and river slope were
   mirrored north–south against `ibound`, `dem`, `landmask`, `min_dem`, `Z0` and
   the capillary parameters. This package aligns every field geographically,
   using a verified per-file table (`STORAGE_ORDER` in `src/inputs.jl`).
   `audit_static_inputs` re-checks the orientation on every build. Numerical
   agreement with the accepted reference run is therefore **not expected**,
   although the record counts still match (1,410 CHD, 15,165 RIV).
2. **MODFLOW 6 version.** 6.3.0, the version Modflow6_jll provides. The
   version of the reference `mf6` executable is not recorded.
3. **Record sets.** RCH covers every active cell (a zero rate adds no flow). The
   Python code skipped zero-rate cells. EVT is restricted to active cells, which
   excludes no cells for this dataset.

The Python reference behaviours below are kept on purpose: the IMS settings
(`outer_dvclose = 5 m`), DRN disabled, NaN → 0 in restart heads, the global
steady head as the first IC and CHD head, and scaling RCH/EVT by
true area / DIS area.

## Not yet done

- Human case: capped pumping (`pumpingCapacity_...nc`, reference Case 3). WEL
  and its validation exist; the capacity capping does not.
- The VIC side: routing forcing and the soil-moisture update from
  `capillary_rise` (`vic_runner.update_statefile`) belong in the Julia VIC.
