#=
Port of the parameter transformations in plugin/config_module.py. NaN semantics
follow numpy: `max`/`min` propagate NaN (np.maximum/np.minimum) and comparisons
with NaN are false (np.where). Line references are to config_module.py.
=#

const RHO_WATER = 1000.0
const MU_WATER = 0.001
const G_GRAVITY = 9.81
const RIV_MANNING = 0.045
const RIV_RESISTANCE = 1.0
const DEFAULT_AQDEPTH = 200.0
const SECONDS_PER_DAY = 86400.0

"Intrinsic permeability (log10 m²) to hydraulic conductivity (m/day)."
permeability_to_conductivity(log10k) = RHO_WATER * 10.0^log10k * (G_GRAVITY / MU_WATER) * 86400.0

"np.clip(x, a_min=lo) keeping NaN."
clip_min(x, lo) = isnan(x) ? x : max(x, lo)

struct RiverParameters
    cells::Vector{CartesianIndex{2}}
    is_river::BitVector             # riv_cond > 0 (width > 30 m); otherwise drain-like
    width::Vector{Float64}          # m
    slope::Vector{Float64}          # usable slope, NaN when ≤ 5e-5
    bottom_bankfull::Vector{Float64}
    head_ini::Vector{Float64}       # fallback stage
    conductance::Vector{Float64}    # m²/day
    bottom::Vector{Float64}         # MF6 RBOT
end

struct StaticParameters
    top::Matrix{Float64}
    botm::Array{Float64,3}
    k::Array{Float64,3}
    k33::Array{Float64,3}
    ss::Array{Float64,3}
    initial_head::Array{Float64,3}
    chd_cells::Vector{CartesianIndex{2}}
    chd_head::Matrix{Float64}             # (ncell, nlay)
    rch_cells::Vector{CartesianIndex{2}}
    area_factor::Matrix{Float64}          # true cell area / DIS cell area
    river::RiverParameters
    evt_cells::Vector{CartesianIndex{2}}
    evt_surface::Vector{Float64}
    evt_rate::Vector{Float64}             # area-corrected, m/day
    evt_depth::Vector{Float64}
    wel_cells::Vector{CartesianIndex{2}}
end

"Layer geometry (`cal_aqdepth`, `cal_toplayer_elevation`, `cal_botlayer_elevation`)."
function layer_geometry(s::StaticInputs)
    dem = [lm ? d : 0.0 for (lm, d) in zip(s.landmask, s.dem)]
    top = [isnan(d) ? 0.0 : d for d in dem]
    aqdepth = map(zip(s.landmask, s.aqdepth)) do (lm, a)
        a = lm ? a : DEFAULT_AQDEPTH
        (isnan(a) || !(a >= 0)) ? DEFAULT_AQDEPTH : a
    end
    botm = cat(top .- 0.1 .* aqdepth, top .- aqdepth; dims = 3)
    return top, aqdepth, botm
end

"NPF conductivities and STO storage coefficient (`get_npf_param`)."
function hydraulic_parameters(s::StaticInputs, aqdepth)
    ratio = s.cell_area ./ (s.delr * s.delc)
    thick1 = [lm ? 0.1a : NaN for (lm, a) in zip(s.landmask, aqdepth)]
    thick2 = [lm ? a - 0.1a : NaN for (lm, a) in zip(s.landmask, aqdepth)]
    k_inp = permeability_to_conductivity.(s.ksat_log)
    kh1_ini = clip_min.(ifelse.(isnan.(s.ksat_l1_conf_log), k_inp, permeability_to_conductivity.(s.ksat_l1_conf_log)), 0.01)
    kh2_ini = clip_min.(ifelse.(isnan.(s.ksat_l2_conf_log), k_inp, permeability_to_conductivity.(s.ksat_l2_conf_log)), 0.01)

    kv1_ini = ifelse.(kh1_ini .> -999, kh1_ini .* ratio, 10.0 .* ratio)
    kv2 = 10.0 .* ratio
    kv1 = max.(thick1 ./ 5000, kv1_ini)
    kh1 = max.(30.0, kh1_ini .* thick1) ./ thick1      # transmissivity ≥ 30 m²/day
    kh2 = max.(30.0, kh2_ini .* thick2) ./ thick2
    kh1[kh1 .< 0] .= 20.0
    kh2[kh2 .< 0] .= 190.0
    kv1[kv1 .< 0] .= 1e10
    kv2[kv2 .< 0] .= 1e10
    nan_to(a, v) = [isnan(x) ? v : x for x in a]
    k = cat(nan_to(kh1, 20.0), nan_to(kh2, 20.0); dims = 3)
    k33 = cat(nan_to(kv1, 20.0), nan_to(kv2, 20.0); dims = 3)

    # `missingvalue` in the reference is aqdepth_ini[0][0], which is NaN for Indus.
    sy = [lm ? x : NaN for (lm, x) in zip(s.landmask, s.spe_yi)]
    sy = max.(min.(max.(0.01, sy), 1.0), 0.11)
    sy = nan_to(sy, 0.1)
    ss = cat(sy, sy; dims = 3)
    return k, k33, ss
end

"Static river geometry and conductance (`get_riv_param` without the discharge terms)."
function river_parameters(s::StaticInputs, top)
    min_dem2 = [x < 0 ? 0.0 : x for x in s.min_dem]
    z0_1 = max.(min_dem2, s.z0_floodplain)
    z0_2 = ifelse.(s.z0_floodplain .< 0, z0_1, s.z0_floodplain)

    width = 4.8 .* sqrt.(s.qbank)
    slope = ifelse.(width .> 0, s.riv_slope, NaN)
    slope_used = ifelse.(slope .> 0.00005, slope, NaN)

    head_ini = ifelse.(width .> 30, z0_2, top)
    head_ini[.!s.landmask] .= 0.0
    depth_bankfull = ((RIV_MANNING .* sqrt.(s.qbank)) ./ (4.8 .* sqrt.(slope_used))) .^ (3 / 5)
    bottom_bankfull = min_dem2 .- depth_bankfull

    riv_cond = ifelse.(width .> 30, width .* sqrt.(2 .* s.cell_area), 0.0) ./ RIV_RESISTANCE
    riv_cond[.!s.landmask] .= 0.0
    drn_cond = ifelse.(riv_cond .== 0, max.(10.0, width) .* sqrt.(2 .* s.cell_area) ./ RIV_RESISTANCE, 0.0)

    bottom = ifelse.(riv_cond .> 0, bottom_bankfull, head_ini)
    bottom[isinf.(bottom)] .= 0.0
    cond = ifelse.(riv_cond .> 0, riv_cond, drn_cond)

    # Stage is finite on land for any discharge, so the record set is static.
    cells = [c for c in CartesianIndices(s.idomain)
             if s.landmask[c] && s.idomain[c] > 0 && isfinite(bottom[c]) &&
                isfinite(cond[c]) && cond[c] != 0]
    return RiverParameters(
        cells, riv_cond[cells] .> 0, width[cells], slope_used[cells],
        bottom_bankfull[cells], head_ini[cells], cond[cells], bottom[cells],
    )
end

"""
    river_stage(river, i, discharge) -> Float64

River stage of record `i` for a monthly mean discharge in m³/s (dynamic part of
`get_riv_param`). Negative discharge sets the stage to the river bottom so the
river cannot recharge the aquifer.
"""
function river_stage(r::RiverParameters, i::Integer, discharge::Real)
    q = isfinite(discharge) ? Float64(discharge) : 0.0
    q < 0 && return r.bottom[i]
    r.is_river[i] || return r.head_ini[i]
    depth = (RIV_MANNING * q / (r.width[i] * sqrt(r.slope[i])))^(3 / 5)
    depth = isnan(depth) || depth < 0.01 ? 0.0 : depth
    stage = r.bottom_bankfull[i] + depth
    return isnan(stage) ? r.head_ini[i] : stage
end

"EVT records (`get_cpr_param`), restricted to active cells."
function evt_parameters(s::StaticInputs, area_factor)
    cells = [c for c in CartesianIndices(s.idomain)
             if s.landmask[c] && s.idomain[c] > 0 && isfinite(s.cpsurface[c]) &&
                isfinite(s.cprate[c]) && isfinite(s.cpdistinct[c]) && s.cprate[c] != 0]
    return cells, s.cpsurface[cells], s.cprate[cells] .* area_factor[cells], s.cpdistinct[cells]
end

"CHD records in both layers at `idomain == 2` (`get_chd_input`)."
function chd_parameters(s::StaticInputs)
    cells = findall(==(2), s.idomain)
    head = [s.steady_head[c, l] for c in cells, l in 1:2]
    bad = count(!isfinite, head)
    bad == 0 || error("$bad CHD cells have no finite global steady-state head")
    return cells, head
end

"""
    read_initial_head(path, inputs) -> Array{Float64,3}

Restart head `(lon, lat, layer)` from a NetCDF with variable `gwl` or
`groundwater_head` on the ibound grid. NaNs become 0 (Python reference behaviour).
"""
function read_initial_head(path, s::StaticInputs)
    NCDataset(path) do ds
        name = haskey(ds, "gwl") ? "gwl" : "groundwater_head"
        v = ds[name]
        latname = haskey(ds, "lat") ? "lat" : "latitude"
        lat = Float64.(ds[latname][:])
        maximum(abs.(lat .- s.lat)) < 1e-3 || error("$path: latitude order differs from the model grid")
        h = to_float.(v[:, :, :])
        size(h) == (size(s)..., 2) || error("$path:$name has size $(size(h))")
        replace!(h, NaN => 0.0)
    end
end

"""
    build_parameters(inputs, config) -> StaticParameters

All static MODFLOW package data. The first-period initial head is the global
steady-state head with NaN → 0, unless `config.initial_head_file` is set.
"""
function build_parameters(s::StaticInputs, config::Config)
    top, aqdepth, botm = layer_geometry(s)
    k, k33, ss = hydraulic_parameters(s, aqdepth)
    area_factor = s.cell_area ./ (s.delr * s.delc)
    initial_head = if isempty(config.initial_head_file)
        replace(s.steady_head, NaN => 0.0)
    else
        read_initial_head(config.initial_head_file, s)
    end
    chd_cells, chd_head = chd_parameters(s)
    evt = config.foc ? evt_parameters(s, area_factor) :
          (CartesianIndex{2}[], Float64[], Float64[], Float64[])
    return StaticParameters(
        top, botm, k, k33, ss, initial_head,
        chd_cells, chd_head,
        findall(>(0), s.idomain), area_factor,
        river_parameters(s, top),
        evt...,
        config.pumping ? findall(==(1), s.idomain) : CartesianIndex{2}[],
    )
end

"""
    audit_static_inputs(inputs, params) -> Vector{Pair{String,Bool}}

Case 0 static audit (reference_cases.md) plus orientation fingerprints.
"""
function audit_static_inputs(s::StaticInputs, p::StaticParameters)
    active = s.idomain .> 0
    rowmean_area = vec(sum(s.cell_area; dims = 1)) ./ size(s.cell_area, 1)
    land = s.landmask .& isfinite.(s.steady_head[:, :, 1]) .& isfinite.(s.dem)
    sukkur = CartesianIndex(argmin(abs.(s.lon .- 68.87)), argmin(abs.(s.lat .- 27.7)))
    window = CartesianIndices((sukkur[1]-2:sukkur[1]+2, sukkur[2]-2:sukkur[2]+2))
    checks = [
        "grid is 180 rows x 204 columns" => size(s) == (204, 180),
        "idomain only contains 0, 1, 2" => all(in((0, 1, 2)), s.idomain),
        "all idomain == 2 cells have finite heads in both layers" =>
            all(isfinite, s.steady_head[s.idomain .== 2, :]),
        "true cell area positive on active cells" => all(>(0), s.cell_area[active]),
        "DIS top above both layer bottoms on active cells" =>
            all(p.top[active] .> p.botm[:, :, 1][active] .> p.botm[:, :, 2][active]),
        "package arrays finite" => all(isfinite, p.k) && all(isfinite, p.k33) && all(isfinite, p.ss) &&
            all(isfinite, p.botm) && all(isfinite, p.initial_head),
        "orientation: cell area decreases northward" => issorted(rowmean_area; rev = true),
        "orientation: DEM correlates with global steady head" =>
            cor(s.dem[land], s.steady_head[:, :, 1][land]) > 0.8,
        "orientation: bankfull discharge > 5000 m3/s near Sukkur" => maximum(s.qbank[window]) > 5000,
        "orientation: bankfull discharge > 10 m3/s only on land" => all(s.landmask[s.qbank .> 10]),
    ]
    return checks
end
