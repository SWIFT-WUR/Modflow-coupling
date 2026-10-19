#=
Grid convention
---------------
All 2D arrays are `(ncol, nrow) = (lon, lat)` with latitude ascending: the grid
of the `ibound` input, sorted south-first. Julia's column-major order then
equals MODFLOW 6 DIS node order (column fastest), so MODFLOW row 1 is the
southernmost row. With a uniform DIS without georeference the mirrored y-axis
has no physical effect.

Every other input is selected at the ibound cell centres by coordinate, so it
may cover a larger area (a basin cut from global maps). Latitude labels are
trusted unless the config sets `lat_order` for that input: some datasets carry
wrong labels (examples/indus/README.md), and `audit_static_inputs` checks the
orientation with physical fingerprints on every model build.
=#
"Static inputs aligned to the model grid. Missing values are NaN."
struct StaticInputs
    lon::Vector{Float64}
    lat::Vector{Float64}
    delr::Float64
    delc::Float64
    idomain::Matrix{Int}
    landmask::BitMatrix
    dem::Matrix{Float64}
    aqdepth::Matrix{Float64}
    cell_area::Matrix{Float64}
    ksat_log::Matrix{Float64}
    ksat_l1_conf_log::Matrix{Float64}
    ksat_l2_conf_log::Matrix{Float64}
    spe_yi::Matrix{Float64}
    qbank::Matrix{Float64}
    riv_slope::Matrix{Float64}
    min_dem::Matrix{Float64}
    z0_floodplain::Matrix{Float64}
    cpsurface::Matrix{Float64}
    cprate::Matrix{Float64}
    cpdistinct::Matrix{Float64}
    kq3::Matrix{Float64}            # groundwater recession constant, 1/day
    catchment::Matrix{Float64}      # catchment id
    steady_head::Array{Float64,3}   # NaN without a `steady_head` input
end

Base.size(s::StaticInputs) = size(s.idomain)

"Convert raw NetCDF values to Float64, mapping missing, non-finite and fill values (≥ 1e30) to NaN."
to_float(x) = ismissing(x) || !isfinite(x) || abs(x) >= 1e30 ? NaN : Float64(x)

const LON_NAMES = ("lon", "longitude", "x")
const LAT_NAMES = ("lat", "latitude", "y")

"Name of the data variable of `spec`: the given one, or the only non-coordinate variable with ≥ 2 dims."
function data_variable(ds, spec::InputSpec, path)
    isempty(spec.variable) || return spec.variable
    names = [k for k in keys(ds) if ndims(ds[k]) >= 2 && !(k in LON_NAMES) && !(k in LAT_NAMES)]
    length(names) == 1 || error("$path holds $(length(names)) data variables ($(join(names, ", "))): set `variable`")
    return only(names)
end

"""
File coordinates of variable `v`: longitudes, and latitudes in storage order
(sorted per `lat_order` when the labels are untrusted).
"""
function file_coordinates(ds, v, path, lat_order::Symbol)
    dims = dimnames(v)
    (length(dims) >= 2 && dims[1] in LON_NAMES && dims[2] in LAT_NAMES) ||
        error("$path:$(name(v)) has dims $dims, expected (lon, lat, ...) in Julia order")
    flon = Float64.(ds[dims[1]][:])
    flat = Float64.(ds[dims[2]][:])
    lat_order === :south_first && sort!(flat)
    lat_order === :north_first && sort!(flat; rev = true)
    return flon, flat
end

"""
    read_field(path, spec, lon, lat; tol = 1e-3) -> Array{Float64}

Read the variable of `spec` (dims `(lon, lat[, layer])`) at the grid
coordinates `lon`, `lat` by nearest-coordinate selection. Every grid coordinate
must have a file coordinate within `tol` degrees; matching shapes alone are not
accepted as proof of identical cell order.
"""
function read_field(path, spec::InputSpec, lon, lat; tol = 1e-3)
    NCDataset(path) do ds
        v = ds[data_variable(ds, spec, path)]
        flon, flat = file_coordinates(ds, v, path, spec.lat_order)
        ii = [nearest_index(flon, x; tol) for x in lon]
        jj = [nearest_index(flat, y; tol) for y in lat]
        (i1, i2), (j1, j2) = extrema(ii), extrema(jj)
        rest = ntuple(_ -> Colon(), ndims(v) - 2)
        block = v[i1:i2, j1:j2, rest...]
        nlay = size(block, 3)
        a = [to_float(block[i - i1 + 1, j - j1 + 1, l]) for i in ii, j in jj, l in 1:nlay]
        ndims(v) == 2 ? a[:, :, 1] : a
    end
end

"Index of the nearest coordinate; ties go to the larger index (pandas `nearest`)."
function nearest_index(coords::AbstractVector, x; tol)
    best, dist = 0, Inf
    for (i, c) in pairs(coords)
        d = abs(c - x)
        if d <= dist
            best, dist = i, d
        end
    end
    dist <= tol || error("no coordinate within $tol of $x (nearest is $dist away)")
    return best
end

"""
Model grid from the ibound input: longitudes (ascending) and latitudes sorted
south-first, and the uniform DIS spacing in metres (degrees × 111 km, as
`config.delrow/delcol` of the Python reference).
"""
function read_grid(path, spec::InputSpec)
    NCDataset(path) do ds
        v = ds[data_variable(ds, spec, path)]
        lon, lat = file_coordinates(ds, v, path, spec.lat_order)
        issorted(lon) || error("$path: longitudes must be ascending")
        lat = sort(lat)
        dx, dy = (lon[end] - lon[1]) / (length(lon) - 1), (lat[end] - lat[1]) / (length(lat) - 1)
        lon, lat, dx * 111_000.0, dy * 111_000.0
    end
end

"""
    read_static_inputs(config) -> StaticInputs

Read all static inputs of `config.inputs` and align them to the grid of the
`ibound` input. Optional inputs that are not given are NaN.
"""
function read_static_inputs(config::Config)
    validate(config)
    lon, lat, delr, delc = read_grid(input_path(config, "ibound"), config.inputs["ibound"])
    nan_field = fill(NaN, length(lon), length(lat))
    field(key) = haskey(config.inputs, key) ?
        read_field(input_path(config, key), config.inputs[key], lon, lat) : copy(nan_field)
    idomain = field("ibound")
    all(isfinite, idomain) || error("ibound has missing values on the grid")
    # The steady head may come from a coarser-labelled global file (e.g. Float32
    # coordinates): accept the nearest cell centre within 0.6 grid cells.
    steady_head = haskey(config.inputs, "steady_head") ?
        read_field(input_path(config, "steady_head"), config.inputs["steady_head"], lon, lat;
                   tol = 0.6 * min(lon[2] - lon[1], lat[2] - lat[1])) :
        fill(NaN, length(lon), length(lat), 2)
    size(steady_head, 3) == 2 || error("steady_head has $(size(steady_head, 3)) layers, expected 2")

    return StaticInputs(
        lon, lat, delr, delc, Int.(idomain),
        map(x -> isfinite(x) && x != 0, field("landmask")),
        field("dem"),
        field("aqdepth"),
        field("cell_area"),
        field("ksat_log"),
        field("ksat_l1_conf_log"),
        field("ksat_l2_conf_log"),
        field("spe_yi"),
        field("qbank"),
        field("riv_slope"),
        field("min_dem"),
        field("z0_floodplain"),
        field("cpsurface"),
        field("cprate"),
        field("cpdistinct"),
        field("kq3"),
        field("catchment"),
        steady_head,
    )
end

"""
    align_to_grid(field, lon, lat, inputs; atol = 1e-3) -> Matrix{Float64}

Reorder a `(lon, lat)` field (e.g. VIC output) to the model grid by
nearest-coordinate selection, as `mfrun._align_vic_field` does. Every model
coordinate must have a source coordinate within `atol` degrees; matching
shapes alone are not accepted as proof of identical cell order.
"""
function align_to_grid(field::AbstractMatrix, lon::AbstractVector, lat::AbstractVector,
                       inputs::StaticInputs; atol = 1e-3)
    size(field) == (length(lon), length(lat)) ||
        error("field size $(size(field)) does not match coordinates $((length(lon), length(lat)))")
    ii = [nearest_index(lon, x; tol = atol) for x in inputs.lon]
    jj = [nearest_index(lat, y; tol = atol) for y in inputs.lat]
    return [to_float(field[i, j]) for i in ii, j in jj]
end
