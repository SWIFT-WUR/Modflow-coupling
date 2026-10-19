"Path of the forcing file for the month starting at `date`."
forcing_path(config::Config, date::Date) =
    replace(config.forcing_file, "{yyyy}" => Dates.format(date, "yyyy"), "{mm}" => Dates.format(date, "mm"))

"""
    read_forcing(config, inputs, date) -> (recharge, discharge)

Read the month containing `date` from the land-surface model output (VIC
`OUT_GWRECHARGE` and `OUT_DISCHARGE` by default, see `Config`) and convert it as
`vic_runner.PostProcessVIC` does: recharge mm/month → m/day, discharge stays
m3/s. Fields are aligned to the model grid by coordinate. A file with a single time step is used as is;
otherwise the time step within the month is selected.
"""
function read_forcing(config::Config, s::StaticInputs, date::Date)
    path = forcing_path(config, date)
    NCDataset(path) do ds
        coord(names) = Float64.(ds[first(filter(n -> haskey(ds, n), names))][:])
        lon, lat = coord(LON_NAMES), coord(LAT_NAMES)
        times = haskey(ds, "time") ? ds["time"][:] : [date]
        k = length(times) == 1 ? 1 :
            findfirst(t -> year(t) == year(date) && month(t) == month(date), times)
        k === nothing && error("$path has no time step in $(Dates.format(date, "yyyy-mm"))")
        slice(name) = (v = ds[name]; ndims(v) == 3 ? v[:, :, k] : v[:, :])
        recharge = align_to_grid(slice(config.forcing_recharge_variable), lon, lat, s) ./ daysinmonth(date) ./ 1000
        discharge = align_to_grid(slice(config.forcing_discharge_variable), lon, lat, s)
        recharge, discharge
    end
end
