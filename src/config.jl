"Static input file names inside `data_dir`, keyed by role."
const DEFAULT_FILES = Dict(
    "ibound" => "ibound_indus_normalized.nc",
    "clone" => "Clone_05min_indus.nc",
    "dem" => "dem_ini_indus.nc",
    "aqdepth" => "aqdepth_ini_indus.nc",
    "cell_area" => "cell_area_indus.nc",
    "ksat_log" => "ksat_log_indus.nc",
    "ksat_l1_conf_log" => "ksat_l1_conf_log_indus.nc",
    "ksat_l2_conf_log" => "ksat_l2_conf_log_indus.nc",
    "spe_yi" => "spe_yi_inp_indus.nc",
    "landmask" => "landmask_indus.nc",
    "qbank" => "bankfull_discharge_indus.nc",
    "riv_slope" => "riv_slope_indus.nc",
    "min_dem" => "min_dem_indus.nc",
    "z0_floodplain" => "Z0_floodplain_indus_updated.nc",
    "capillary" => "capillary_rise_indus.nc",
    "kq3" => "KQ3_indus.nc",
    "ldd" => "ldd_indus_updated.nc",
    "global_steady_head" => "mf_ss_gwl_5minDRN.nc",
)

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

- `foc`: fully online coupling; adds the EVT package and capillary-rise feedback.
- `pumping`: adds a WEL package on all `idomain == 1` cells (human extension).
- `forcing_file`: optional VIC monthly output (OUT_GWRECHARGE in mm/month,
  OUT_DISCHARGE in m3/s). `{yyyy}` and `{mm}` in the path are replaced by the
  stress-period year and month. When set, each `update` reads the month's
  fields, unless they were set through BMI since the previous update.
- `initial_head_file`: optional restart head (variable `gwl` or `groundwater_head`,
  dims `(layer, lat, lon)` on the ibound grid). NaNs are replaced by 0 as in the
  Python reference. When empty, the global steady-state head is used.
"""
Base.@kwdef struct Config
    data_dir::String
    work_dir::String
    output_dir::String
    start_date::Date
    end_date::Date
    period::String = "month"
    model_name::String = "INDUS"
    foc::Bool = true
    pumping::Bool = false
    initial_head_file::String = ""
    forcing_file::String = ""
    write_monthly_output::Bool = true
    save_mf6_binary::Bool = false
    files::Dict{String,String} = copy(DEFAULT_FILES)
    solver::SolverSettings = SolverSettings()
end

"""
    load_config(path) -> Config

Read a TOML run configuration. Relative paths are resolved against the
directory of the TOML file. See `scripts/indus_natural_foc.toml` for an example.
"""
function load_config(path::AbstractString)
    toml = TOML.parsefile(path)
    root = dirname(abspath(path))
    resolve(p) = isempty(p) ? p : normpath(isabspath(p) ? p : joinpath(root, p))

    paths = toml["paths"]
    time = toml["time"]
    coupling = get(toml, "coupling", Dict{String,Any}())
    output = get(toml, "output", Dict{String,Any}())
    files = merge(DEFAULT_FILES, Dict{String,String}(get(toml, "files", Dict{String,Any}())))
    solver_table = get(toml, "solver", Dict{String,Any}())
    solver = SolverSettings(; (Symbol(k) => v for (k, v) in solver_table)...)

    config = Config(;
        data_dir = resolve(paths["data_dir"]),
        work_dir = resolve(paths["work_dir"]),
        output_dir = resolve(paths["output_dir"]),
        start_date = Date(time["start_date"]),
        end_date = Date(time["end_date"]),
        period = lowercase(get(time, "period", "month")),
        model_name = uppercase(get(toml, "model_name", "INDUS")),
        foc = get(coupling, "foc", true),
        pumping = get(coupling, "pumping", false),
        initial_head_file = resolve(get(paths, "initial_head_file", "")),
        forcing_file = resolve(get(paths, "forcing_file", "")),
        write_monthly_output = get(output, "write_monthly", true),
        save_mf6_binary = get(output, "save_mf6_binary", false),
        files,
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

input_path(config::Config, key::AbstractString) = joinpath(config.data_dir, config.files[key])
