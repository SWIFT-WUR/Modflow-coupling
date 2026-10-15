"""MODFLOW -> mGV coupling example.

How to run:
    julia --project=examples/mgv_coupling examples/mgv_coupling/couple_mgv.jl \
        <mgv_config.toml> [modflow_config.toml]

This reproduces the two MODFLOW -> VIC feedbacks of the Python plugin:

  baseflow_to_routing [m3/s]     -> mGV `channel_water__area~cross-sectional`
      Python: `discharge_mf`, read by VIC-WUR's routing plugin as MONTH forcing.
      Here: held for the month; before every daily mGV step one day of baseflow
      volume is added to the channel storage, which mGV's kinematic wave then
      releases downstream.
  capillary_rise      [mm/month] -> mGV `soil__moisture_layer_N`, bottom layer (FOC only)
      Python: vic_runner.update_statefile, added to STATE_SOIL_MOISTURE layer 3
      once at the start of the month. Here: spread evenly over the days of the
      month and added before every daily mGV step. Excess above
      `soil__moisture~maximum` moves up one layer at a time.

Timing: mGV runs month M with the feedback from MODFLOW month M-1, then
MODFLOW runs month M, then its feedback is set for mGV month M+1.

The VIC -> MODFLOW direction is NOT coupled here: MODFLOW reads its own offline
`forcing_file` (see examples/synthetic). mGV has no groundwater-recharge output
yet, so there is nothing equivalent to OUT_GWRECHARGE to send.
"""
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
Add a capillary-rise depth [mm] to the bottom soil layer, then move any excess
above `soil__moisture~maximum` up one layer at a time. As in the Python plugin,
only cells that received capillary rise are capped, and water that does not fit
in the top layer is dropped. Returns that dropped volume [mm summed over cells].
"""
function add_capillary_rise!(vic, cpr::Vector{Float32})
    n = length(cpr)
    nlayers = length(BMI.get_value_ptr(vic, "soil__moisture")) ÷ n
    # one column per layer, the bottom layer last
    w = reduce(hcat, (get_field(vic, "soil__moisture_layer_$k", n) for k in 1:nlayers))
    wmax = reduce(hcat, (get_field(vic, "soil__moisture~maximum_layer_$k", n) for k in 1:nlayers))
    # mGV holds NaN outside its land domain; leave those cells alone
    active = (cpr .> 0) .& vec(all(isfinite, w; dims = 2)) .& vec(all(isfinite, wmax; dims = 2))
    w[active, nlayers] .+= cpr[active]
    overflow = 0.0
    for k in nlayers:-1:1
        excess = @. ifelse(active, max(w[:, k] - wmax[:, k], 0f0), 0f0)
        w[:, k] .-= excess
        k > 1 ? (w[:, k-1] .+= excess) : (overflow = sum(Float64, excess))
    end
    for k in 1:nlayers
        BMI.set_value(vic, "soil__moisture_layer_$k", w[:, k])
    end
    return overflow
end

"""
Add one mGV time step of baseflow [m3/s] to the channel storage. mGV stores the
water in a channel as its wetted cross-section `channel_water__area~cross-sectional`
[m2], so a volume `baseflow * dt` raises it by `baseflow * dt / channel__length`. The
area cannot go below zero, so a loss (negative baseflow) larger than the water
in the channel is only partly taken. Returns the volume applied [m3].
"""
function add_baseflow!(vic, baseflow::Vector{Float32}, cell_length::Vector{Float32}, dt)
    area = get_field(vic, "channel_water__area~cross-sectional", length(baseflow))
    new_area = @. ifelse(cell_length > 0, max(area + baseflow * dt / cell_length, 0f0), area)
    BMI.set_value(vic, "channel_water__area~cross-sectional", new_area)
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
length(BMI.get_value_ptr(vic, "channel_water__area~cross-sectional")) == n_vic ||
    error("mGV routing grid differs from its land grid")
dt = Float32(BMI.get_time_step(vic))   # s
# NaN where the routing input has no channel; those cells take no baseflow
cell_length = replace(get_field(vic, "channel__length", n_vic), NaN32 => 0f0)

mf_to_vic = cell_map(mf_x, mf_y, vic_x, vic_y)
n_shared = count(>(0), mf_to_vic)
@info "Grid overlap" mGV_cells = n_vic MODFLOW_cells = length(mf_x) * length(mf_y) shared = n_shared
n_shared == 0 && @warn "mGV and MODFLOW grids do not overlap: the feedback will be zero"

baseflow = zeros(Float32, n_vic)   # m3/s on the mGV grid; MODFLOW feedback of the previous month

# --- monthly loop ------------------------------------------------------------
while (
    BMI.get_current_time(gw) < BMI.get_end_time(gw)
    && lastdayofmonth(current_date(gw)) <= mgv_end
)
    month = current_date(gw)
    
    cpr = regrid(BMI.get_value_ptr(gw, "capillary_rise"), mf_to_vic) / daysinmonth(month)
    
    for _ in 1:daysinmonth(month)
        add_capillary_rise!(vic, cpr)
        add_baseflow!(vic, baseflow, cell_length, dt)
        BMI.update(vic)
    end

    BMI.update(gw)

    baseflow .= regrid(BMI.get_value_ptr(gw, "baseflow_to_routing"), mf_to_vic)

end

BMI.finalize(gw)
BMI.finalize(vic)
