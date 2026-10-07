#=
MODFLOW -> mGV coupling example: IndusModflow feeds mGV through BMI, in memory.

    julia --project=examples/mgv_coupling examples/mgv_coupling/couple_mgv.jl \
        <mgv_config.toml> [modflow_config.toml]

This reproduces the two MODFLOW -> VIC feedbacks of the Python plugin
(coupling_contract.md):

  baseflow_to_routing [m3/s]     -> mGV `river_channel_area` (channel storage)
      Python: `discharge_mf`, read by VIC-WUR's routing plugin as MONTH forcing.
      Here: held for the month; before every daily mGV step one day of baseflow
      volume is added to the channel storage, which mGV's kinematic wave then
      releases downstream.
  capillary_rise      [mm/month] -> mGV `moisture`, bottom layer (FOC only)
      Python: vic_runner.update_statefile, added to STATE_SOIL_MOISTURE layer 3
      once at the start of the month. Excess above `maximum_moisture` moves up
      one layer at a time.

Timing: mGV runs month M with the feedback from MODFLOW month M-1, then
MODFLOW runs month M, then its feedback is set for mGV month M+1.

The VIC -> MODFLOW direction is NOT coupled here: MODFLOW reads its own offline
`forcing_file` (see examples/synthetic). mGV has no groundwater-recharge output
yet, so there is nothing equivalent to OUT_GWRECHARGE to send.

mGV needs the BMI branch `11-requested-feature-implement-bmi` plus
mgv_bmi_coupling_inputs.patch (adds `river_channel_area`, `river_cell_length`,
`river_discharge` and `maximum_moisture` to its BMI).
=#
using IndusModflow, Dates, Printf
import BasicModelInterface as BMI
import mGV

const COORD_ATOL = 1e-3    # degrees; both models are on the same 5 arcmin grid

mgv_config = abspath(get(ARGS, 1, ""))
mf_config = abspath(get(ARGS, 2, joinpath(@__DIR__, "..", "synthetic", "synthetic.toml")))
isfile(mgv_config) || error("usage: couple_mgv.jl <mgv_config.toml> [modflow_config.toml]")

"""
For every destination cell, the linear index of the source cell with the same
centre (within `atol`), or 0 when the source grid does not cover it. Both
grids are flattened longitude-fastest, as both BMIs return them.
"""
function cell_map(src_x, src_y, dst_x, dst_y; atol = COORD_ATOL)
    near(v, x) = (i = argmin(abs.(v .- x)); abs(v[i] - x) <= atol ? i : 0)
    ix = [near(src_x, x) for x in dst_x]
    iy = [near(src_y, y) for y in dst_y]
    nx = length(src_x)
    return [(ix[i] > 0 && iy[j] > 0) ? ix[i] + (iy[j] - 1) * nx : 0
            for i in eachindex(dst_x), j in eachindex(dst_y)] |> vec
end

"Source values on the destination grid; uncovered and non-finite cells become 0."
function regrid(src, map)
    return Float32[(k > 0 && isfinite(src[k])) ? src[k] : 0f0 for k in map]
end

get_field(model, name, n) = BMI.get_value(model, name, zeros(Float32, n))

"""
Add a monthly capillary-rise depth [mm] to the bottom soil layer, then move any
excess above `maximum_moisture` up one layer at a time. As in the Python plugin,
only cells that received capillary rise are capped, and water that does not fit
in the top layer is dropped. Returns that dropped volume [mm summed over cells].
"""
function add_capillary_rise!(vic, cpr::Vector{Float32})
    n = length(cpr)
    nlayers = length(BMI.get_value_ptr(vic, "moisture")) ÷ n
    # (lon, lat, layer) flattened: one column per layer, the bottom layer last
    w = reshape(get_field(vic, "moisture", n * nlayers), n, nlayers)
    wmax = reshape(get_field(vic, "maximum_moisture", n * nlayers), n, nlayers)
    # mGV holds NaN outside its land domain; leave those cells alone
    active = (cpr .> 0) .& vec(all(isfinite, w; dims = 2)) .& vec(all(isfinite, wmax; dims = 2))
    w[active, nlayers] .+= cpr[active]
    overflow = 0.0
    for k in nlayers:-1:1
        excess = @. ifelse(active, max(w[:, k] - wmax[:, k], 0f0), 0f0)
        w[:, k] .-= excess
        k > 1 ? (w[:, k-1] .+= excess) : (overflow = sum(Float64, excess))
    end
    BMI.set_value(vic, "moisture", vec(w))
    return overflow
end

"""
Add one mGV time step of baseflow [m3/s] to the channel storage. mGV stores the
water in a channel as its wetted cross-section `river_channel_area` [m2], so a
volume `baseflow * dt` raises it by `baseflow * dt / river_cell_length`. The
area cannot go below zero, so a loss (negative baseflow) larger than the water
in the channel is only partly taken. Returns the volume applied [m3].
"""
function add_baseflow!(vic, baseflow::Vector{Float32}, cell_length::Vector{Float32}, dt)
    area = get_field(vic, "river_channel_area", length(baseflow))
    new_area = @. ifelse(cell_length > 0, max(area + baseflow * dt / cell_length, 0f0), area)
    BMI.set_value(vic, "river_channel_area", new_area)
    return sum(Float64, @. ifelse(cell_length > 0, (new_area - area) * cell_length, 0f0))
end

# --- initialize --------------------------------------------------------------
gw = BMI.initialize(IndusGroundwater, mf_config)
vic = BMI.initialize(mGV.Model, mgv_config)

BMI.get_time_step(vic) == 86400 || error("the example steps mGV per day")
mgv_start = Date(unix2datetime(BMI.get_start_time(vic)))
mgv_end = Date(unix2datetime(BMI.get_end_time(vic)))
mgv_start == current_date(gw) ||
    error("mGV starts $mgv_start, MODFLOW starts $(current_date(gw))")

mf_x = BMI.get_grid_x(gw, 0, zeros(BMI.get_grid_shape(gw, 0)[2]))
mf_y = BMI.get_grid_y(gw, 0, zeros(BMI.get_grid_shape(gw, 0)[1]))
vic_x = BMI.get_grid_x(vic, 0)
vic_y = BMI.get_grid_y(vic, 0)
n_vic = length(vic_x) * length(vic_y)
length(BMI.get_value_ptr(vic, "river_channel_area")) == n_vic ||
    error("mGV routing grid differs from its land grid")
dt = Float32(BMI.get_time_step(vic))   # s
# NaN where the routing input has no channel; those cells take no baseflow
cell_length = replace(get_field(vic, "river_cell_length", n_vic), NaN32 => 0f0)

mf_to_vic = cell_map(mf_x, mf_y, vic_x, vic_y)
n_shared = count(>(0), mf_to_vic)
@info "Grid overlap" mGV_cells = n_vic MODFLOW_cells = length(mf_x) * length(mf_y) shared = n_shared
n_shared == 0 && @warn "mGV and MODFLOW grids do not overlap: the feedback will be zero"

baseflow = zeros(Float32, n_vic)   # m3/s on the mGV grid; MODFLOW feedback of the previous month

# --- monthly loop ------------------------------------------------------------
while BMI.get_current_time(gw) < BMI.get_end_time(gw) &&
      lastdayofmonth(current_date(gw)) <= mgv_end
    month = current_date(gw)

    # 1. mGV month M, one update per day. (mGV's update_until does not work yet
    #    on this branch.) MODFLOW's baseflow from M-1 goes into the channel
    #    storage before every step; capillary rise was added once, below.
    applied = 0.0
    for _ in 1:daysinmonth(month)
        applied += add_baseflow!(vic, baseflow, cell_length, dt)
        BMI.update(vic)
    end
    requested = sum(Float64, @. ifelse(cell_length > 0, baseflow, 0f0)) * dt * daysinmonth(month)
    Date(unix2datetime(BMI.get_current_time(vic))) == lastdayofmonth(month) ||
        error("mGV is out of step with MODFLOW at $month")

    # 2. MODFLOW month M (inputs from its offline forcing_file, see header)
    BMI.update(gw)

    # 3. MODFLOW -> mGV, for month M+1
    #    Positive baseflow drains the aquifer into the river; negative takes river
    #    water into the aquifer. mGV still routes its own subsurface_runoff too:
    #    in a two-way coupling, groundwater must not reach the river twice.
    baseflow .= regrid(BMI.get_value_ptr(gw, "baseflow_to_routing"), mf_to_vic)

    #    capillary_rise is zero without FOC, which leaves the soil unchanged.
    cpr = regrid(BMI.get_value_ptr(gw, "capillary_rise"), mf_to_vic)
    overflow = add_capillary_rise!(vic, cpr)

    @printf("%s  baseflow into channels %.3g of %.3g m3  |  capillary rise for next month %.1f mm·cells (%.1f did not fit)\n",
            Dates.format(month, "yyyy-mm"), applied, requested, sum(cpr), overflow)
end

BMI.finalize(gw)   # mGV has no BMI.finalize yet; its writer closes at the end of each year
