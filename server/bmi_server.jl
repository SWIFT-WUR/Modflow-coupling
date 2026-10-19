#=
Serve the ModflowCoupling MODFLOW model over Remote BMI (https://github.com/eWaterCycle/remotebmi),
e.g. for the BMI Controller (https://github.com/BSchilperoort/bmi-controller).

    julia --project=server server/bmi_server.jl          # port 50051, or BMI_PORT

Initialize with an absolute config path, e.g.
    /path/to/ModflowCoupling/examples/indus/synthetic.toml
=#
using ModflowCoupling
import RemoteBMI.Server: run_bmi_server

port = parse(Int, get(ENV, "BMI_PORT", "50051"))
@info "ModflowCoupling Remote BMI server on port $port"
run_bmi_server(GroundwaterModel, "0.0.0.0", port)
