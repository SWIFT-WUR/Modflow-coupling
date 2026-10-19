"""
Static inputs, keyed by role as in the `[inputs]` table of a config. Every
field is a `(lon, lat)` map read from a NetCDF file, except `steady_head`
(`(lon, lat, layer)`). `required` is `true`, `false` or `:foc` (needed only
with `foc = true`).
"""
const INPUT_ROLES = [
    # role                 required  description
    ("ibound",             true,     "model domain: 0 inactive, 1 active, 2 constant head (CHD); defines the grid"),
    ("landmask",           true,     "land cells: finite and non-zero"),
    ("dem",                true,     "land surface elevation (m), DIS top"),
    ("aqdepth",            true,     "aquifer thickness (m)"),
    ("cell_area",          true,     "true cell area (m2)"),
    ("ksat_log",           true,     "log10 intrinsic permeability (m2)"),
    ("ksat_l1_conf_log",   false,    "log10 permeability of layer 1 where confined; ksat_log elsewhere"),
    ("ksat_l2_conf_log",   false,    "log10 permeability of layer 2 where confined; ksat_log elsewhere"),
    ("spe_yi",             true,     "specific yield (-)"),
    ("qbank",              true,     "bankfull discharge (m3/s)"),
    ("riv_slope",          true,     "river slope (-)"),
    ("min_dem",            true,     "minimum elevation in the cell (m)"),
    ("z0_floodplain",      true,     "floodplain elevation (m): drain level and river-bottom reference"),
    ("kq3",                true,     "groundwater recession constant (1/day)"),
    ("catchment",          true,     "catchment ids (DRN levels are not below a catchment's lowest Z0)"),
    ("cpsurface",          :foc,     "EVT surface elevation (m)"),
    ("cprate",             :foc,     "EVT maximum rate (m/day)"),
    ("cpdistinct",         :foc,     "EVT extinction depth (m)"),
    ("steady_head",        false,    "two-layer head (m), (lon, lat, layer): initial head and CHD head; may cover a larger grid"),
]
const INPUT_KEYS = first.(INPUT_ROLES)

"""
    InputSpec(file; variable = "", lat_order = :coordinates)

One static input. `file` is relative to `data_dir` (or absolute). With an empty
`variable`, the file must hold exactly one data variable. `lat_order` says how
latitude is stored:

- `:coordinates` (default): trust the file's latitude labels;
- `:south_first` or `:north_first`: the labels are wrong and the data are stored
  in this order (the label values are only used as a set).

Fields are selected at the grid coordinates of `ibound`, so the file may cover a
larger area (e.g. the globe) at the same resolution.
"""
Base.@kwdef struct InputSpec
    file::String
    variable::String = ""
    lat_order::Symbol = :coordinates
end
InputSpec(file::AbstractString; kw...) = InputSpec(; file = String(file), kw...)

const LAT_ORDERS = (:coordinates, :south_first, :north_first)

function InputSpec(key::AbstractString, entry)
    entry isa AbstractString && return InputSpec(entry)
    entry isa AbstractDict || error("inputs.$key must be a file name or a table")
    unknown = setdiff(keys(entry), ("file", "variable", "lat_order"))
    isempty(unknown) || error("inputs.$key: unknown keys $(join(unknown, ", "))")
    haskey(entry, "file") || error("inputs.$key has no `file`")
    order = Symbol(get(entry, "lat_order", "coordinates"))
    order in LAT_ORDERS || error("inputs.$key: lat_order must be one of $(join(LAT_ORDERS, ", "))")
    return InputSpec(entry["file"]; variable = get(entry, "variable", ""), lat_order = order)
end

"IMS settings; defaults are those of the Python reference (`mf_run.mfrun`)."
Base.@kwdef struct SolverSettings
    outer_dvclose::Float64 = 5.0
    inner_dvclose::Float64 = 5.0
    inner_rclose::Float64 = 10.0
    outer_maximum::Int = 200
    inner_maximum::Int = 500
    under_relaxation_gamma::Float64 = 0.99
    relaxation_factor::Float64 = 0.0
end

"""
    Config

Run configuration. Stress periods are calendar months from the month of
`start_date` through the month of `end_date` (inclusive), or with
`period = "day"` single days from `start_date` through `end_date`. The Python
reference uses months. `period = "steady"` is one steady-state stress period:
one `update` solves for the heads in balance with the inputs set through BMI
(rates per day, as for daily periods); `end_date` is then ignored.

- `inputs`: the static input files by role (`INPUT_ROLES`).
- `foc`: fully online coupling; adds the EVT package and capillary-rise feedback.
- `pumping`: adds a WEL package on all `idomain == 1` cells (human extension).
- `forcing_file`: optional monthly land-surface model output (recharge in
  mm/month, discharge in m3/s; variable names `forcing_recharge_variable` and
  `forcing_discharge_variable`, VIC names by default). `{yyyy}` and `{mm}` in
  the path are replaced by the stress-period year and month. When set, each
  `update` reads the month's fields, unless they were set through BMI since the
  previous update.
- `initial_head_file`: optional restart head (variable `gwl` or `groundwater_head`,
  dims `(layer, lat, lon)` on the ibound grid). NaNs are replaced by 0 as in the
  Python reference. When empty, the `steady_head` input is used.
- `check_orientation`: run the orientation fingerprints of `audit_static_inputs`.

`Config(config; field = value, ...)` copies a config with some fields replaced.
"""
Base.@kwdef struct Config
    data_dir::String
    work_dir::String
    output_dir::String
    start_date::Date
    end_date::Date
    inputs::Dict{String,InputSpec}
    period::String = "month"
    model_name::String = "GWF"
    foc::Bool = true
    pumping::Bool = false
    initial_head_file::String = ""
    forcing_file::String = ""
    forcing_recharge_variable::String = "OUT_GWRECHARGE"
    forcing_discharge_variable::String = "OUT_DISCHARGE"
    write_monthly_output::Bool = true
    save_mf6_binary::Bool = false
    check_orientation::Bool = true
    solver::SolverSettings = SolverSettings()
end

Config(c::Config; kw...) = Config(; (f => get(kw, f, getfield(c, f)) for f in fieldnames(Config))...)

const CONFIG_KEYS = ("model_name", "paths", "time", "coupling", "output", "audit", "forcing", "solver", "inputs")

"""
    load_config(path) -> Config

Read a TOML run configuration. One file describes the whole MODFLOW model.
Relative paths in `[paths]` are resolved against the directory of the TOML
file, and input files in `[inputs]` against `data_dir`, so a config and its
data can be moved together. See `examples/indus/indus.toml`.
"""
function load_config(path::AbstractString)
    toml = TOML.parsefile(path)
    unknown = setdiff(keys(toml), CONFIG_KEYS)
    isempty(unknown) || error("$path: unknown keys $(join(sort(collect(unknown)), ", "))")
    root = dirname(abspath(path))
    resolve(p) = isempty(p) ? p : normpath(isabspath(p) ? p : joinpath(root, p))
    haskey(toml, "paths") || error("$path has no [paths] table")
    toml["paths"] = Dict{String,Any}(k => resolve(v) for (k, v) in toml["paths"])
    section(name) = get(toml, name, Dict{String,Any}())
    paths, time, coupling, output, forcing, audit = section.(("paths", "time", "coupling", "output", "forcing", "audit"))
    haskey(toml, "inputs") || error("$path has no [inputs] table")
    inputs = Dict{String,InputSpec}(k => InputSpec(k, v) for (k, v) in toml["inputs"])
    solver = SolverSettings(; (Symbol(k) => v for (k, v) in section("solver"))...)

    config = Config(;
        data_dir = paths["data_dir"],
        work_dir = paths["work_dir"],
        output_dir = paths["output_dir"],
        start_date = Date(time["start_date"]),
        end_date = Date(time["end_date"]),
        inputs,
        period = lowercase(get(time, "period", "month")),
        model_name = uppercase(get(toml, "model_name", "GWF")),
        foc = get(coupling, "foc", true),
        pumping = get(coupling, "pumping", false),
        initial_head_file = get(paths, "initial_head_file", ""),
        forcing_file = get(paths, "forcing_file", ""),
        forcing_recharge_variable = get(forcing, "recharge_variable", "OUT_GWRECHARGE"),
        forcing_discharge_variable = get(forcing, "discharge_variable", "OUT_DISCHARGE"),
        write_monthly_output = get(output, "write_monthly", true),
        save_mf6_binary = get(output, "save_mf6_binary", false),
        check_orientation = get(audit, "orientation", true),
        solver,
    )
    validate(config)
    return config
end

function validate(config::Config)
    config.period in ("month", "day", "steady") ||
        error("period must be \"month\", \"day\" or \"steady\", got \"$(config.period)\"")
    config.period == "month" && day(config.start_date) != 1 &&
        error("start_date must be the first day of a month, got $(config.start_date)")
    config.end_date >= config.start_date || error("end_date is before start_date")
    length(config.model_name) <= 16 || error("MODFLOW 6 model names are limited to 16 characters")
    config.pumping && !config.foc && @warn "pumping without FOC: capillary rise is not returned"

    unknown = setdiff(keys(config.inputs), INPUT_KEYS)
    isempty(unknown) || error("unknown inputs $(join(sort(collect(unknown)), ", ")); known roles: $(join(INPUT_KEYS, ", "))")
    missing_inputs = [key for (key, required, _) in INPUT_ROLES
                      if (required === true || (required === :foc && config.foc)) && !haskey(config.inputs, key)]
    isempty(missing_inputs) || error("missing inputs: $(join(missing_inputs, ", "))")
    isempty(config.initial_head_file) && !haskey(config.inputs, "steady_head") &&
        error("set paths.initial_head_file or inputs.steady_head for the initial head")
    return config
end

"Length of one stress period (a steady-state period counts as one day)."
period_step(config::Config) = config.period == "month" ? Month(1) : Day(1)

"First day of every stress period."
stress_periods(config::Config) =
    config.period == "steady" ? [config.start_date] :
    config.period == "day" ? collect(config.start_date:Day(1):config.end_date) :
    collect(config.start_date:Month(1):firstdayofmonth(config.end_date))

"Number of days in the stress period starting at `date`."
period_days(config::Config, date::Date) = Dates.value(date + period_step(config) - date)

input_path(config::Config, key::AbstractString) = joinpath(config.data_dir, config.inputs[key].file)
