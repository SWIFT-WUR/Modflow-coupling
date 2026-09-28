#=
BMI for the coupled exchange (coupling_contract.md). Grid 0 is the 2D model grid
(lat × lon, uniform rectilinear); grid 1 adds the two layers. Values are
flattened in Julia order: longitude fastest, latitude ascending, then layer.
Values sit at cell centres (location "node"; grid origin is the SW cell centre).
Time is in days since `start_date`, reported as CF units "days since <start_date>";
one `update` is one calendar month.
=#

const INPUT_VARS = Dict(
    "groundwater_recharge" => "m d-1",
    "river_discharge" => "m3 s-1",
    "groundwater_abstraction" => "m3 d-1",
)
const OUTPUT_VARS = Dict(
    "groundwater_head" => "m",
    "groundwater_head_layer_1" => "m",
    "groundwater_head_layer_2" => "m",
    "groundwater_depth" => "m",
    "baseflow_to_routing" => "m3 s-1",
    "capillary_rise" => "mm month-1",
)

function BMI.initialize(::Type{IndusGroundwater}, config_file::AbstractString)
    return IndusGroundwater(load_config(config_file))
end

function BMI.update(m::IndusGroundwater)
    advance!(m)
    return nothing
end

function BMI.update_until(m::IndusGroundwater, time::Float64)
    while BMI.get_current_time(m) < time - 1e-9
        advance!(m)
    end
    return nothing
end

function BMI.finalize(m::IndusGroundwater)
    finalize!(m)
    return nothing
end

BMI.get_component_name(::IndusGroundwater) = "IndusModflow (MODFLOW 6)"
BMI.get_input_item_count(::IndusGroundwater) = length(INPUT_VARS)
BMI.get_output_item_count(::IndusGroundwater) = length(OUTPUT_VARS)
BMI.get_input_var_names(::IndusGroundwater) = sort!(collect(keys(INPUT_VARS)))
BMI.get_output_var_names(::IndusGroundwater) = sort!(collect(keys(OUTPUT_VARS)))

function var_units(name)
    haskey(INPUT_VARS, name) && return INPUT_VARS[name]
    haskey(OUTPUT_VARS, name) && return OUTPUT_VARS[name]
    error("unknown variable $name")
end

BMI.get_var_units(::IndusGroundwater, name::String) = var_units(name)
BMI.get_var_grid(::IndusGroundwater, name::String) = (var_units(name); name == "groundwater_head" ? 1 : 0)
BMI.get_var_type(::IndusGroundwater, name::String) = (var_units(name); "Float64")
BMI.get_var_itemsize(::IndusGroundwater, name::String) = (var_units(name); sizeof(Float64))
BMI.get_var_nbytes(m::IndusGroundwater, name::String) = length(var_array(m, name)) * sizeof(Float64)
BMI.get_var_location(::IndusGroundwater, name::String) = (var_units(name); "node")

BMI.get_start_time(::IndusGroundwater) = 0.0
BMI.get_end_time(m::IndusGroundwater) = Float64(Dates.value(m.periods[end] + Month(1) - m.periods[1]))
BMI.get_current_time(m::IndusGroundwater) = Float64(Dates.value(current_date(m) - m.periods[1]))
BMI.get_time_units(m::IndusGroundwater) = "days since $(m.periods[1])"
"Length of the next stress period (days)."
BMI.get_time_step(m::IndusGroundwater) = Float64(daysinmonth(current_date(m)))

function var_array(m::IndusGroundwater, name::String)
    name == "groundwater_recharge" && return m.recharge
    name == "river_discharge" && return m.discharge
    name == "groundwater_abstraction" && return m.abstraction
    name == "groundwater_head" && return m.head
    name == "groundwater_head_layer_1" && return view(m.head, :, :, 1)
    name == "groundwater_head_layer_2" && return view(m.head, :, :, 2)
    # Depth of the layer-1 head below land surface (derived; a copy, not a pointer)
    name == "groundwater_depth" && return ifelse.(m.inputs.landmask, m.params.top .- m.head[:, :, 1], NaN)
    name == "baseflow_to_routing" && return m.baseflow
    name == "capillary_rise" && return m.capillary
    error("unknown variable $name")
end

BMI.get_value_ptr(m::IndusGroundwater, name::String) = vec(var_array(m, name))
BMI.get_value(m::IndusGroundwater, name::String, dest) = copyto!(dest, vec(var_array(m, name)))
BMI.get_value_at_indices(m::IndusGroundwater, name::String, dest, inds) =
    copyto!(dest, vec(var_array(m, name))[inds])

function BMI.set_value(m::IndusGroundwater, name::String, value)
    length(value) == length(m.recharge) || error("$name needs $(length(m.recharge)) values, got $(length(value))")
    name == "groundwater_recharge" ? set_recharge!(m, value) :
    name == "river_discharge" ? set_discharge!(m, value) :
    name == "groundwater_abstraction" ? set_abstraction!(m, value) :
    error("$name is not an input variable")
    return nothing
end

function BMI.set_value_at_indices(m::IndusGroundwater, name::String, inds, value)
    full = copy(vec(var_array(m, name)))
    full[inds] .= value
    BMI.set_value(m, name, full)
end

BMI.get_grid_rank(::IndusGroundwater, grid::Int) = grid == 0 ? 2 : 3
BMI.get_grid_size(m::IndusGroundwater, grid::Int) = length(m.recharge) * (grid == 0 ? 1 : 2)
BMI.get_grid_type(::IndusGroundwater, ::Int) = "uniform_rectilinear"

function BMI.get_grid_shape(m::IndusGroundwater, grid::Int)
    ncol, nrow = size(m.inputs)
    return grid == 0 ? [nrow, ncol] : [2, nrow, ncol]
end

BMI.get_grid_shape(m::IndusGroundwater, grid::Int, shape::AbstractVector{Int}) =
    copyto!(shape, BMI.get_grid_shape(m, grid))

function BMI.get_grid_spacing(m::IndusGroundwater, grid::Int, spacing::AbstractVector{Float64})
    s = m.inputs
    dlat, dlon = s.lat[2] - s.lat[1], s.lon[2] - s.lon[1]
    return copyto!(spacing, grid == 0 ? [dlat, dlon] : [1.0, dlat, dlon])
end

function BMI.get_grid_origin(m::IndusGroundwater, grid::Int, origin::AbstractVector{Float64})
    s = m.inputs
    return copyto!(origin, grid == 0 ? [s.lat[1], s.lon[1]] : [0.0, s.lat[1], s.lon[1]])
end

BMI.get_grid_x(m::IndusGroundwater, ::Int, x::AbstractVector{Float64}) = copyto!(x, m.inputs.lon)
BMI.get_grid_y(m::IndusGroundwater, ::Int, y::AbstractVector{Float64}) = copyto!(y, m.inputs.lat)
