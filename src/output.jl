const TIME_UNITS = "days since 1979-01-01 00:00:00"

"""
    write_monthly_output(model)

Append the completed month to `output_dir/indus_modflow_monthly.nc` and write a
restart head `output_dir/head_YYYYMMDD.nc` (end-of-month date, variable `gwl`,
Float64) that `initial_head_file` accepts.
"""
function write_monthly_output(m::IndusGroundwater)
    dir = m.config.output_dir
    mkpath(dir)
    s = m.inputs
    start = m.periods[m.completed]
    stop = lastdayofmonth(start)

    path = joinpath(dir, "indus_modflow_monthly.nc")
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

    NCDataset(joinpath(dir, "head_$(Dates.format(stop, "yyyymmdd")).nc"), "c") do ds
        defVar(ds, "lon", s.lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", s.lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "layer", [0, 1], ("layer",))
        defVar(ds, "gwl", m.head, ("lon", "lat", "layer");
               attrib = ["units" => "m", "long_name" => "groundwater head at end of $stop"],
               fillvalue = NaN)
    end
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
           attrib = ["units" => "m3 s-1", "long_name" => "-(RIV flow)/86400, VIC routing forcing for the next month"], fillvalue = NaN)
    defVar(ds, "capillary_rise", Float64, grid2;
           attrib = ["units" => "mm month-1", "long_name" => "EVT outflow depth, added to the next VIC month's state"], fillvalue = NaN)
    defVar(ds, "groundwater_recharge", Float64, grid2;
           attrib = ["units" => "m d-1", "long_name" => "VIC recharge input for the stress period"], fillvalue = NaN)
    defVar(ds, "river_discharge", Float64, grid2;
           attrib = ["units" => "m3 s-1", "long_name" => "VIC discharge input for the stress period"], fillvalue = NaN)
    return ds
end
