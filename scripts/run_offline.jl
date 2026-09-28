#=
Run IndusModflow offline from monthly VIC output named by `paths.forcing_file`
in the config (e.g. the reference run's VIC output, or the synthetic forcing
from scripts/make_synthetic_forcing.jl).

    julia --project=. scripts/run_offline.jl examples/synthetic/synthetic.toml
=#
using IndusModflow, Dates
import BasicModelInterface as BMI

function main(config_file)
    model = BMI.initialize(IndusGroundwater, config_file)
    isempty(model.config.forcing_file) && error("$config_file sets no paths.forcing_file")
    try
        while BMI.get_current_time(model) < BMI.get_end_time(model)
            date = current_date(model)
            BMI.update(model)
            println(Dates.format(date, "yyyy-mm"), ": baseflow_to_routing sum ",
                    round(sum(filter(isfinite, model.baseflow)), digits = 1), " m3/s")
        end
    finally
        BMI.finalize(model)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1])
end
