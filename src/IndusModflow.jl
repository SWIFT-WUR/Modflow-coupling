"""
    IndusModflow

MODFLOW 6 groundwater component for the Indus VIC–MODFLOW coupling, driven
in-memory through ModflowInterface.jl. Replaces the file-based Python plugin
(`plugin/`); formulas follow `indus_vic_julia_bmi_modflow_handoff/coupling_contract.md`.
"""
module IndusModflow

using Dates
using Statistics: cor
using TOML
using NCDatasets
using ModflowInterface
import ModflowInterface as MF
import BasicModelInterface as BMI

export IndusGroundwater, Config, load_config, read_static_inputs, build_parameters,
       audit_static_inputs, align_to_grid, current_date, package_flows

include("config.jl")
include("inputs.jl")
include("parameters.jl")
include("mf6_writer.jl")
include("forcing.jl")
include("model.jl")
include("output.jl")
include("bmi.jl")

end
