"""
    ModflowCoupling

MODFLOW 6 groundwater component for coupling with a land-surface model (mGV,
VIC) through BMI, driven in-memory through ModflowInterface.jl. The domain and
its static inputs come from the `[inputs]` table of a TOML config; the domain
can be a single basin or the globe. Replaces the file-based Python plugin
(`plugin/`); formulas follow `coupling_contract.md` of the Indus VIC–MODFLOW handoff.
"""
module ModflowCoupling

using Dates
using Statistics: cor, mean
using TOML
using NCDatasets
using ModflowInterface
import ModflowInterface as MF
import BasicModelInterface as BMI

export GroundwaterModel, Config, InputSpec, load_config, read_static_inputs, build_parameters,
       audit_static_inputs, align_to_grid, current_date, package_flows, write_head

include("config.jl")
include("inputs.jl")
include("parameters.jl")
include("mf6_writer.jl")
include("forcing.jl")
include("model.jl")
include("output.jl")
include("bmi.jl")

end
