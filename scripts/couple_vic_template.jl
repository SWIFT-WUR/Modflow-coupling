#=
Template for the monthly VIC–MODFLOW loop with a BMI-enabled Julia VIC.
The VIC-side names (`VICModel`, its variable names, month stepping) are
placeholders: replace them once the Julia VIC BMI is fixed.

Timing (coupling_contract.md): VIC runs month M, MODFLOW runs month M with
VIC's month-M fields, and MODFLOW's month-M feedback is consumed by VIC in M+1.
=#
using IndusModflow, Dates
import BasicModelInterface as BMI

# using VIC                                           # placeholder
# vic = BMI.initialize(VIC.VICModel, "vic_indus_natural_foc.toml")
gw = BMI.initialize(IndusGroundwater, joinpath(@__DIR__, "indus_natural_foc.toml"))

# VIC grid coordinates, used to align exchanged fields by coordinate (not by shape)
# vic_lon = BMI.get_grid_x(vic, grid, zeros(nlon)); vic_lat = BMI.get_grid_y(vic, grid, zeros(nlat))

while BMI.get_current_time(gw) < BMI.get_end_time(gw)
    month = current_date(gw)

    # 1. VIC month M (it already holds the feedback from M-1, see step 4)
    # BMI.update_until(vic, <end of month in VIC time units>)

    # 2. VIC -> MODFLOW: monthly means on the MODFLOW grid
    # recharge  = align_to_grid(<OUT_GWRECHARGE mm/month> ./ daysinmonth(month) ./ 1000, vic_lon, vic_lat, gw.inputs)  # m/day
    # discharge = align_to_grid(<OUT_DISCHARGE m3/s>, vic_lon, vic_lat, gw.inputs)
    # BMI.set_value(gw, "groundwater_recharge", vec(recharge))
    # BMI.set_value(gw, "river_discharge", vec(discharge))

    # 3. MODFLOW month M (one stress period of daysinmonth(month) days)
    BMI.update(gw)

    # 4. MODFLOW -> VIC, for month M+1
    baseflow = BMI.get_value_ptr(gw, "baseflow_to_routing")   # m3/s, NaN outside landmask
    capillary = BMI.get_value_ptr(gw, "capillary_rise")       # mm/month, add once to VIC soil layer 3
    # BMI.set_value(vic, "<routing forcing: discharge_mf>", <baseflow on the VIC grid>)
    # BMI.set_value(vic, "<soil moisture increment>", <capillary on the VIC grid>)  # see vic_runner.update_statefile
end

BMI.finalize(gw)
# BMI.finalize(vic)
