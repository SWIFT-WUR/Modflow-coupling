#=
Write synthetic monthly VIC output for offline Indus runs, in the layout
of VIC image-driver output (lat north-first, lon, time; NaN off the landmask):

- OUT_GWRECHARGE (mm/month): a monsoon peak in Jul–Sep (strongest in the wetter
  north-east), a smaller western-disturbance peak in Feb–Mar, and a 1 mm/month base.
- OUT_DISCHARGE (m3/s): bankfull discharge × a seasonal fraction peaking in
  Jul–Aug (snowmelt + monsoon), 0.12 in winter.
Each year gets a wetness factor so the years differ (1980 wet, 1981 dry).

The fields are illustrative only; they are not a VIC simulation.

    julia --project=. examples/indus/make_synthetic_forcing.jl [output.nc] [first_year] [last_year]
=#
using ModflowCoupling, NCDatasets, Dates

gauss(x, mu, sigma) = exp(-0.5 * ((x - mu) / sigma)^2)

function make_synthetic_forcing(path, years; config = load_config(joinpath(@__DIR__, "indus.toml")))
    s = read_static_inputs(config)
    lon, lat = s.lon, s.lat
    wetness(y) = (0.9, 1.3, 0.7)[mod(y - first(years), 3) + 1]
    # Spatial monsoon weight: rises towards the north-east (Punjab, Himalayan foothills)
    monsoon_weight = [clamp(0.15 + 0.5 * (la - 23) / 15 + 0.6 * (lo - 66) / 17, 0.1, 1.2) for lo in lon, la in lat]

    dates = [Date(y, m, 15) for y in years for m in 1:12]
    recharge = fill(NaN32, length(lon), length(lat), length(dates))
    discharge = fill(NaN32, length(lon), length(lat), length(dates))
    for (k, d) in enumerate(dates)
        m, w = month(d), wetness(year(d))
        seasonal_q = 0.12 + 0.88 * gauss(m, 7.5, 1.5)
        for c in CartesianIndices(s.landmask)
            s.landmask[c] || continue
            r = 1.0 + w * (60.0 * monsoon_weight[c] * gauss(m, 8.0, 1.0) + 15.0 * gauss(m, 2.5, 0.8))
            recharge[c, k] = r
            discharge[c, k] = s.qbank[c] * seasonal_q * w
        end
    end

    mkpath(dirname(path))
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        # VIC image output stores latitude north-first; ModflowCoupling aligns by coordinate.
        defVar(ds, "lon", lon, ("lon",); attrib = ["units" => "degrees_east", "standard_name" => "longitude"])
        defVar(ds, "lat", reverse(lat), ("lat",); attrib = ["units" => "degrees_north", "standard_name" => "latitude"])
        defVar(ds, "time", DateTime.(dates), ("time",); attrib = ["units" => "days since 1979-01-01 00:00:00", "calendar" => "standard"])
        defVar(ds, "OUT_GWRECHARGE", reverse(recharge; dims = 2), ("lon", "lat", "time");
               attrib = ["units" => "mm/month", "long_name" => "synthetic groundwater recharge"],
               fillvalue = NaN32, deflatelevel = 4)
        defVar(ds, "OUT_DISCHARGE", reverse(discharge; dims = 2), ("lon", "lat", "time");
               attrib = ["units" => "m3/s", "long_name" => "synthetic river discharge"],
               fillvalue = NaN32, deflatelevel = 4)
        ds.attrib["title"] = "Synthetic VIC monthly output for offline Indus MODFLOW runs"
        ds.attrib["history"] = "Created $(now()) by ModflowCoupling/examples/indus/make_synthetic_forcing.jl"
    end
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = get(ARGS, 1, joinpath(@__DIR__, "vic_synthetic_monthly.nc"))
    y0 = parse(Int, get(ARGS, 2, "1979"))
    y1 = parse(Int, get(ARGS, 3, "1981"))
    println("wrote ", make_synthetic_forcing(normpath(out), y0:y1))
end
