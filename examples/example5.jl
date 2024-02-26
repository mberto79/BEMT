using DelimitedFiles
using Printf

d_string = @sprintf "%.2f" 15e-2 
base = "examples/data/example3/"

# Example providing path with joinpath
RPM_file = readdlm(joinpath(base,"RPM_d$(d_string).csv"))
T_file = readdlm(joinpath(base,"T_d$(d_string).csv"))
Q_file = readdlm(joinpath(base,"Q_d$(d_string).csv"))


# Plot results
p1 = plot(
    RPM_file, T_file, 
    label="D = $d_string m", xlabel="RPM", ylabel="T [N]"
    )
p2 = plot(
    RPM_file, Q_file,
    label="D = $d_string m", xlabel="RPM", ylabel="Q [Nm]"
    )


d_string = @sprintf "%.2f" 25e-2 
base = "examples/data/example4/"

# Example providing path with joinpath
RPM_file = readdlm(joinpath(base,"RPM_d$(d_string).csv"))
T_file = readdlm(joinpath(base,"T_d$(d_string).csv"))
Q_file = readdlm(joinpath(base,"Q_d$(d_string).csv"))


# Plot results
plot!(p1, RPM_file, T_file, label="D = $d_string m")
plot!(p2, RPM_file, Q_file, label="D = $d_string m")

# Plot with some plot attributes (see Plots.jl docs)
plot(
    p1,p2, 
    framestyle=:box,
    fg_legend=:false,
    layout=(2,1),
    size=(600,500)
    ) 

savefig("examples/example4_combining_results.svg")