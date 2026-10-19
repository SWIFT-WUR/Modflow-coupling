#= Shared functions of the mGV <-> MODFLOW coupling, used by couple_mgv.jl and
spinup_steady_state.jl. See couple_mgv.jl for what is exchanged and when. =#
using ModflowCoupling, Dates, Printf
import BasicModelInterface as BMI
import mGV

const COORD_ATOL = 1e-3    # degrees; both models are on the same 5 arcmin grid

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

"""
Cell maps between the MODFLOW grid (`gw`) and the mGV grid (`vic`), by cell
centre. Returns `(mf_to_vic, vic_to_mf, n_vic)`: for every mGV cell the MODFLOW
cell index, for every MODFLOW cell the mGV cell index (0 where not covered, see
`cell_map`), and the number of mGV cells. Warns about cells without a match.
"""
function grid_maps(gw, vic)
    mf_x = BMI.get_grid_x(gw, 0, zeros(BMI.get_grid_shape(gw, 0)[2]))
    mf_y = BMI.get_grid_y(gw, 0, zeros(BMI.get_grid_shape(gw, 0)[1]))
    vic_x = BMI.get_grid_x(vic, 0)
    vic_y = BMI.get_grid_y(vic, 0)
    n_vic = length(vic_x) * length(vic_y)
    length(BMI.get_value_ptr(vic, "channel_water__area~cross-sectional")) == n_vic ||
        error("mGV routing grid differs from its land grid")

    mf_to_vic = cell_map(mf_x, mf_y, vic_x, vic_y)
    n_shared = count(>(0), mf_to_vic)
    @info "Grid overlap" mGV_cells = n_vic MODFLOW_cells = length(mf_x) * length(mf_y) shared = n_shared
    n_shared == 0 && @warn "mGV and MODFLOW grids do not overlap: the feedback will be zero"
    vic_to_mf = cell_map(vic_x, vic_y, mf_x, mf_y)
    n_missing = count(==(0), vic_to_mf)
    n_missing > 0 && @warn "$n_missing MODFLOW cells are outside the mGV grid: they get zero recharge and discharge"
    return mf_to_vic, vic_to_mf, n_vic
end

"""
State of a daily mGV <-> MODFLOW exchange: the grid maps, the channel lengths
and MODFLOW's feedback from the previous day. `requested` and `applied` are
the baseflow volumes [m3] over all days so far; their difference is the river
loss that mGV's channels could not supply (water created by the coupling).
"""
mutable struct Coupling
    mf_to_vic::Vector{Int}
    vic_to_mf::Vector{Int}
    n_vic::Int
    dt::Float32                    # s
    cell_length::Vector{Float32}   # m; 0 where mGV has no channel
    baseflow::Vector{Float32}      # m3/s on the mGV grid, from the previous MODFLOW day
    cpr::Vector{Float32}           # mm/day on the mGV grid, from the previous MODFLOW day
    requested::Float64
    applied::Float64
end

"Check that `gw` and `vic` can be coupled daily, and set up the exchange."
function Coupling(gw, vic)
    BMI.get_time_step(vic) == 86400 || error("the coupling steps mGV per day")
    BMI.get_time_step(gw) == 1 || error("MODFLOW must run daily periods: set `period = \"day\"` under [time]")
    mgv_start = Date(unix2datetime(BMI.get_start_time(vic)))
    mgv_start == current_date(gw) ||
        error("mGV starts $mgv_start, MODFLOW starts $(current_date(gw))")
    # Without this, mGV also routes its subsurface runoff to the river, and the
    # recharge sent to MODFLOW would reach the river twice.
    vic.config.API.external_groundwater_routing ||
        error("set `external_groundwater_routing = true` in the [API] section of the mGV config")
    mf_to_vic, vic_to_mf, n_vic = grid_maps(gw, vic)
    # NaN where the routing input has no channel; those cells take no baseflow
    cell_length = replace(get_field(vic, "channel__length", n_vic), NaN32 => 0f0)
    return Coupling(mf_to_vic, vic_to_mf, n_vic, Float32(BMI.get_time_step(vic)), cell_length,
                    zeros(Float32, n_vic), zeros(Float32, n_vic), 0.0, 0.0)
end

"Last day both models can run."
last_day(gw, vic) = min(Date(unix2datetime(BMI.get_end_time(vic))),
                        current_date(gw) + Day(round(Int, BMI.get_end_time(gw) - BMI.get_current_time(gw))) - Day(1))

"""
Run one coupled day D: mGV with MODFLOW's feedback from day D-1, then MODFLOW
with mGV's values of day D. Returns `(discharge, recharge)` as sent to MODFLOW,
on the MODFLOW grid (m3/s and m/day; MODFLOW applies only positive recharge).
"""
function step_day!(c::Coupling, gw, vic)
    # 1. mGV day D, with MODFLOW's feedback from day D-1
    add_capillary_rise!(vic, c.cpr)
    c.requested += sum(Float64, @. ifelse(c.cell_length > 0, c.baseflow, 0f0)) * c.dt
    c.applied += add_baseflow!(vic, c.baseflow, c.cell_length, c.dt)
    BMI.update(vic)

    # 2. mGV day D -> MODFLOW day D
    discharge = Float64.(regrid(get_field(vic, "channel_water__volume_flow_rate", c.n_vic), c.vic_to_mf))  # m3/s
    recharge = Float64.(regrid(get_field(vic, "land_subsurface__runoff", c.n_vic) ./ 1000, c.vic_to_mf))  # mm/day -> m/day
    BMI.set_value(gw, "river_discharge", discharge)
    BMI.set_value(gw, "groundwater_recharge", recharge)
    BMI.update(gw)

    # 3. MODFLOW day D -> mGV day D+1
    c.baseflow .= regrid(BMI.get_value_ptr(gw, "baseflow_to_routing"), c.mf_to_vic)
    c.cpr .= regrid(BMI.get_value_ptr(gw, "capillary_rise"), c.mf_to_vic)
    return discharge, recharge
end

"River loss [m3] that mGV's channels could not supply so far."
created_water(c::Coupling) = c.applied - c.requested
