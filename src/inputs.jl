#=
Grid convention
---------------
All 2D arrays are `(ncol, nrow) = (lon, lat)` with latitude ascending, i.e. the
storage order of `ibound_indus_normalized.nc`. Julia's column-major order then
equals MODFLOW 6 DIS node order (column fastest), so MODFLOW row 1 is the
southernmost row. With a uniform DIS without georeference the mirrored y-axis
has no physical effect.

Latitude storage order of the static inputs
-------------------------------------------
The static NetCDF files do NOT share one latitude storage order, and two of them
(`dem_ini_indus.nc`, `landmask_indus.nc`) carry descending latitude labels while
their data are stored south-first. Coordinate labels can therefore not be used to
align them. The table below records the verified storage order of each field
(2026-09-28), established with physical fingerprints independent of the labels:

- cell area must decrease northward (cos(lat));
- elevation fields must correlate with the geographically referenced global
  steady-state head (r ≈ 0.96 south-first, negative when mirrored);
- `qbank` must peak on the lower Indus main stem (Sukkur, Hyderabad);
- `riv_slope` must correlate with the DEM gradient (r = 0.79 vs -0.26);
- sea cells (SW corner) must be missing in `aqdepth`, `ksat_log`, `spe_yi`.

`audit_static_inputs` repeats the cheap fingerprints on every model build.
The Python reference combined raw `.values` and therefore mirrored the
north-first fields against `ibound`.
=#
const STORAGE_ORDER = Dict(
    "ibound" => :south_first,
    "dem" => :south_first,          # lat labels say north-first: labels are wrong
    "landmask" => :south_first,     # lat labels say north-first: labels are wrong
    "min_dem" => :south_first,
    "z0_floodplain" => :south_first,
    "capillary" => :south_first,    # cpsurface, cprate, cpdistinct (max_moist is north-first)
    "aqdepth" => :north_first,
    "cell_area" => :north_first,
    "ksat_log" => :north_first,
    "ksat_l1_conf_log" => :north_first,
    "ksat_l2_conf_log" => :north_first,
    "spe_yi" => :north_first,
    "qbank" => :north_first,
    "riv_slope" => :north_first,
)

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
    steady_head::Array{Float64,3}
end

Base.size(s::StaticInputs) = size(s.idomain)

"Convert raw NetCDF values to Float64, mapping missing, non-finite and fill values (≥ 1e30) to NaN."
to_float(x) = ismissing(x) || !isfinite(x) || abs(x) >= 1e30 ? NaN : Float64(x)

"""
    read_field(path, var, order, lon, lat) -> Matrix{Float64}

Read a `(lat, lon)` NetCDF variable as a `(lon, lat)` array in model-grid order.
Longitudes must match the grid; latitudes are only checked as a set because
some files carry wrong latitude labels (see `STORAGE_ORDER`).
"""
function read_field(path, var, order::Symbol, lon, lat; atol = 1e-3)
    NCDataset(path) do ds
        v = ds[var]
        dimnames(v) == ("lon", "lat") ||
            error("$path:$var has dims $(dimnames(v)), expected (lon, lat) in Julia order")
        flon = Float64.(ds["lon"][:])
        flat = Float64.(ds["lat"][:])
        (length(flon), length(flat)) == (length(lon), length(lat)) ||
            error("$path:$var has shape $(size(v)), expected $((length(lon), length(lat)))")
        maximum(abs.(flon .- lon)) <= atol || error("$path: longitudes do not match the ibound grid")
        maximum(abs.(sort(flat) .- lat)) <= atol || error("$path: latitudes do not match the ibound grid")
        a = to_float.(v[:, :])
        order === :north_first && reverse!(a; dims = 2)
        a
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
    select_global_head(path, lon, lat) -> Array{Float64,3}

Two-layer global steady-state head at the grid coordinates by nearest-neighbour
selection, as `config_module._get_non_global_steady_head` does. NaNs are kept.
"""
function select_global_head(path, lon, lat)
    NCDataset(path) do ds
        v = ds["gwd"]
        dimnames(v) == ("longitude", "latitude", "layer") ||
            error("$path:gwd has dims $(dimnames(v))")
        glon = Float64.(ds["longitude"][:])
        glat = Float64.(ds["latitude"][:])
        tol_lon = 0.6 * abs(glon[2] - glon[1])
        tol_lat = 0.6 * abs(glat[2] - glat[1])
        ii = [nearest_index(glon, x; tol = tol_lon) for x in lon]
        jj = [nearest_index(glat, y; tol = tol_lat) for y in lat]
        (i1, i2), (j1, j2) = extrema(ii), extrema(jj)
        block = v[i1:i2, j1:j2, :]
        [to_float(block[i - i1 + 1, j - j1 + 1, l]) for i in ii, j in jj, l in 1:size(block, 3)]
    end
end

"Grid spacing in metres from the clone map, as `config.delrow/delcol` (GDAL geotransform × 111 km)."
function clone_spacing(path)
    NCDataset(path) do ds
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        dx = (lon[end] - lon[1]) / (length(lon) - 1)
        dy = abs(lat[end] - lat[1]) / (length(lat) - 1)
        (dx * 111_000.0, dy * 111_000.0)
    end
end

"""
    read_static_inputs(config) -> StaticInputs

Read and align all static inputs to the `ibound_indus_normalized.nc` grid.
"""
function read_static_inputs(config::Config)
    path(key) = input_path(config, key)
    lon, lat, idomain = NCDataset(path("ibound")) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Int.(ds["ibound"][:, :])
    end
    issorted(lat) || error("ibound latitudes must be ascending")
    issorted(lon) || error("ibound longitudes must be ascending")
    field(key, var = "Band1") = read_field(path(key), var, STORAGE_ORDER[key], lon, lat)
    capillary(var) = read_field(path("capillary"), var, STORAGE_ORDER["capillary"], lon, lat)
    delr, delc = clone_spacing(path("clone"))

    return StaticInputs(
        lon, lat, delr, delc, idomain,
        isfinite.(field("landmask")),
        field("dem"),
        field("aqdepth"),
        field("cell_area"),
        field("ksat_log"),
        field("ksat_l1_conf_log"),
        field("ksat_l2_conf_log"),
        field("spe_yi"),
        field("qbank", "qbank"),
        field("riv_slope"),
        field("min_dem"),
        field("z0_floodplain"),
        capillary("cpsurface"),
        capillary("cprate"),
        capillary("cpdistinct"),
        select_global_head(path("global_steady_head"), lon, lat),
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
