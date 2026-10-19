# Indus example

The Indus basin on a 5 arcmin grid (180 lat × 204 lon), the case of the
Indus VIC–MODFLOW reference handoff. The test suite (`test/runtests.jl`) runs
this case.

Files:

- `indus.toml`: the Indus model. Natural FOC run (reference Case 1/2) with
  daily periods 1975–2004 for the mGV coupling (`examples/mgv_coupling`),
  starting from the coupled steady state. Set `period = "month"` and leave out
  `initial_head_file` for the monthly reference setup that the tests run.
- `synthetic.toml`: the same model, monthly 1979–1981 from the global steady
  state, driven by synthetic forcing.
- `make_synthetic_forcing.jl`: writes the synthetic forcing for `synthetic.toml`.

## Input data

The static inputs are not in this repository. The configs expect them in
`../data`, a `data/` directory next to the repository checkout. The tests also
accept `INDUS_DATA=/path/to/data`.

| Role | File |
| --- | --- |
| `ibound` | `ibound_indus_normalized.nc` |
| `dem`, `aqdepth` | `dem_ini_indus.nc`, `aqdepth_ini_indus.nc` |
| `cell_area` | `cell_area_indus.nc` |
| permeability | `ksat_log_indus.nc`, `ksat_l1_conf_log_indus.nc`, `ksat_l2_conf_log_indus.nc` |
| `spe_yi` | `spe_yi_inp_indus.nc` |
| `landmask` | `landmask_indus.nc` |
| rivers | `bankfull_discharge_indus.nc`, `riv_slope_indus.nc`, `min_dem_indus.nc`, `Z0_floodplain_indus_updated.nc` |
| capillary rise (EVT) | `capillary_rise_indus.nc` (`cpsurface`, `cprate`, `cpdistinct`) |
| drainage (DRN) | `KQ3_indus.nc` (recession constant), `ldd_indus_updated.nc` (catchment ids) |
| `steady_head` | `mf_ss_gwl_5minDRN.nc` (global steady state, cut to the grid) |

`Clone_05min_indus.nc`, which the Python reference used for the grid spacing,
is no longer read: the spacing comes from the ibound coordinates. This changes
heads by less than 1e-6 m.

```bash
julia --project=. examples/indus/make_synthetic_forcing.jl   # writes vic_synthetic_monthly.nc here
julia --project=. scripts/run_offline.jl examples/indus/synthetic.toml
```

The synthetic forcing (1979–1981; not a VIC simulation) has a monsoon recharge
peak in Jul–Sep that is strongest in the north-east, and river discharge as a
seasonal fraction of bankfull discharge. The years are normal, wet and dry.

## Latitude storage order

The static NetCDFs do not share one latitude storage order, and two of them
carry wrong labels. The storage order of each field was verified on 2026-09-28
with physical fingerprints that do not depend on the labels:

- cell area must decrease northward (cos(lat));
- elevation fields must correlate with the geographically referenced global
  steady-state head (r ≈ 0.96 south-first, negative when mirrored);
- `qbank` must peak on the lower Indus main stem (Sukkur, Hyderabad);
- `riv_slope` must correlate with the DEM gradient (r = 0.79 vs −0.26);
- sea cells (SW corner) must be missing in `aqdepth`, `ksat_log`, `spe_yi`;
- the Indus main stem (`qbank` > 1000 m³/s) must lie in one `ldd` catchment
  (0.97 north-first vs 0.76 mirrored, 2026-10-16);
- log `KQ3` must correlate with aquifer depth (r = 0.59 vs 0.18) and permeability.

| Stored south-first | Stored north-first |
| --- | --- |
| ibound, **dem**, **landmask**, min_dem, Z0_floodplain, capillary rise (cpsurface, cprate, cpdistinct) | aqdepth, cell_area, ksat_log, ksat_l1/l2_conf_log, spe_yi, qbank, riv_slope, KQ3, ldd |

The labels match the storage order except in `dem_ini_indus.nc` and
`landmask_indus.nc`, whose data are stored south-first under north-first
labels. `indus.toml` therefore sets `lat_order = "south_first"` for those two;
all other inputs are aligned by their labels. `cpdistinct` holds netCDF default
fill values (9.97e36) without a `_FillValue` attribute; values ≥ 1e30 are read
as missing. The generic fingerprints of `audit_static_inputs` run on every
build. The Indus-specific ones (grid size, Sukkur, main-stem catchment, DRN
coverage) are in `indus_fingerprints` in `test/runtests.jl`.

## Differences from the Python reference

1. **Grid orientation (important).** The Python code combined raw `.values`. As
   a result, aquifer depth, cell area, K, storage, bankfull discharge and river
   slope were mirrored north–south against `ibound`, `dem`, `landmask`,
   `min_dem`, `Z0` and the capillary parameters. Here every field is aligned
   geographically (above). Numerical agreement with the accepted reference run
   is therefore **not expected**, although the record counts still match
   (1,410 CHD, 15,165 RIV).
2. **MODFLOW 6 version.** 6.3.0. The version of the reference `mf6` executable
   is not recorded.
3. **DRN is on.** The Python Indus case left out the DRN package of the global
   model (`get_drn_param`). Here it is always on. The global steady-state head
   was computed with this drain. 14,917 of 15,165 active land cells get a drain.
4. **River bottom from `Z0`.** `min_dem` lies a median 8 m above `Z0` in river
   cells, which put rivers above their own drains and let water short-circuit
   from river to aquifer to drain (1,393 km³/year in a steady state with 172
   km³/year of recharge). `NOTES.md` in the data directory lists this and the
   other issues found in the static inputs.
5. **Record sets.** RCH covers every active cell (a zero rate adds no flow). The
   Python code skipped zero-rate cells. EVT is restricted to active cells, which
   excludes no cells for this dataset.
