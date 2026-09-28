#=
Serve IndusModflow over Remote BMI (https://github.com/eWaterCycle/remotebmi),
e.g. for the BMI Controller (https://github.com/BSchilperoort/bmi-controller).

    julia --project=server server/bmi_server.jl          # port 50051, or BMI_PORT

Initialize with an absolute config path, e.g.
    /home/bart/data/modflow/IndusModflow/examples/synthetic/synthetic.toml
=#
using IndusModflow
import RemoteBMI.Server: run_bmi_server

port = parse(Int, get(ENV, "BMI_PORT", "50051"))
@info "IndusModflow Remote BMI server on port $port"
run_bmi_server(IndusGroundwater, "0.0.0.0", port)
