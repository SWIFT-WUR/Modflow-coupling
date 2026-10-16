const TIME_UNITS = "days since 1979-01-01 00:00:00"

"""
    write_monthly_output(model)

Append the completed stress period to `output_dir/indus_modflow_monthly.nc`
(`indus_modflow_daily.nc` with `period = "day"`) and, at the end of each month,
write a restart head `output_dir/head_YYYYMMDD.nc` (end-of-month date, variable
`gwl`, Float64) that `initial_head_file` accepts.
"""
function write_monthly_output(m::IndusGroundwater)
    dir = m.config.output_dir
    mkpath(dir)
    s = m.inputs
    start = m.periods[m.completed]
    stop = start + period_step(m.config) - Day(1)

    path = joinpath(dir, m.config.period == "day" ? "indus_modflow_daily.nc" : "indus_modflow_monthly.nc")
    m.completed == 1 && isfile(path) && rm(path)
    ds = isfile(path) ? NCDataset(path, "a") : create_monthly_file(path, s)
    try
        n = length(ds["time"]) + 1
        ds["time"][n] = Dates.value(start - Date(1979, 1, 1))
        ds["groundwater_head"][:, :, :, n] = m.head
        ds["baseflow_to_routing"][:, :, n] = m.baseflow
        ds["capillary_rise"][:, :, n] = m.capillary
        ds["groundwater_recharge"][:, :, n] = m.recharge
        ds["river_discharge"][:, :, n] = m.discharge
    finally
        close(ds)
    end

    stop == lastdayofmonth(stop) || return
    write_head(m, joinpath(dir, "head_$(Dates.format(stop, "yyyymmdd")).nc");
               long_name = "groundwater head at end of $stop")
end

"""
    write_head(model, path; long_name)

Write the current heads of both layers to `path` (variable `gwl`, Float64, NaN
on inactive cells), in the form `initial_head_file` accepts.
"""
function write_head(m::IndusGroundwater, path::AbstractString;
                    long_name = "groundwater head", attrib = Pair{String,Any}[])
    s = m.inputs
    mkpath(dirname(abspath(path)))
    NCDataset(path, "c"; attrib = attrib) do ds
        defVar(ds, "lon", s.lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", s.lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "layer", [0, 1], ("layer",))
        defVar(ds, "gwl", m.head, ("lon", "lat", "layer");
               attrib = ["units" => "m", "long_name" => long_name], fillvalue = NaN)
    end
    return path
end

function create_monthly_file(path, s::StaticInputs)
    ds = NCDataset(path, "c")
    defDim(ds, "time", Inf)
    defVar(ds, "lon", s.lon, ("lon",); attrib = ["units" => "degrees_east"])
    defVar(ds, "lat", s.lat, ("lat",); attrib = ["units" => "degrees_north"])
    defVar(ds, "layer", [0, 1], ("layer",))
    defVar(ds, "time", Float64, ("time",);
           attrib = ["units" => TIME_UNITS, "calendar" => "standard",
                     "long_name" => "first day of the stress period"])
    grid2 = ("lon", "lat", "time")
    defVar(ds, "groundwater_head", Float64, ("lon", "lat", "layer", "time");
           attrib = ["units" => "m", "long_name" => "head at end of stress period"], fillvalue = NaN)
    defVar(ds, "baseflow_to_routing", Float64, grid2;
           attrib = ["units" => "m3 s-1", "long_name" => "-(RIV flow)/86400, VIC routing forcing for the next period"], fillvalue = NaN)
    defVar(ds, "capillary_rise", Float64, grid2;
           attrib = ["units" => "mm per stress period", "long_name" => "EVT outflow depth, added to the next VIC period's state"], fillvalue = NaN)
    defVar(ds, "groundwater_recharge", Float64, grid2;
           attrib = ["units" => "m d-1", "long_name" => "VIC recharge input for the stress period"], fillvalue = NaN)
    defVar(ds, "river_discharge", Float64, grid2;
           attrib = ["units" => "m3 s-1", "long_name" => "VIC discharge input for the stress period"], fillvalue = NaN)
    return ds
end
