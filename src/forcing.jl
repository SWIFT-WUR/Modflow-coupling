"Monthly VIC output variables read from `forcing_file`."
const VIC_RECHARGE = "OUT_GWRECHARGE"    # mm/month
const VIC_DISCHARGE = "OUT_DISCHARGE"    # m3/s

"Path of the forcing file for the month starting at `date`."
forcing_path(config::Config, date::Date) =
    replace(config.forcing_file, "{yyyy}" => Dates.format(date, "yyyy"), "{mm}" => Dates.format(date, "mm"))

"""
    read_forcing(config, inputs, date) -> (recharge, discharge)

Read the VIC month containing `date` and convert it as `vic_runner.PostProcessVIC`
does: recharge mm/month → m/day, discharge stays m3/s. Fields are aligned to the
model grid by coordinate. A file with a single time step is used as is;
otherwise the time step within the month is selected.
"""
function read_forcing(config::Config, s::StaticInputs, date::Date)
    path = forcing_path(config, date)
    NCDataset(path) do ds
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        times = haskey(ds, "time") ? ds["time"][:] : [date]
        k = length(times) == 1 ? 1 :
            findfirst(t -> year(t) == year(date) && month(t) == month(date), times)
        k === nothing && error("$path has no time step in $(Dates.format(date, "yyyy-mm"))")
        slice(name) = (v = ds[name]; ndims(v) == 3 ? v[:, :, k] : v[:, :])
        recharge = align_to_grid(slice(VIC_RECHARGE), lon, lat, s) ./ daysinmonth(date) ./ 1000
        discharge = align_to_grid(slice(VIC_DISCHARGE), lon, lat, s)
        recharge, discharge
    end
end
