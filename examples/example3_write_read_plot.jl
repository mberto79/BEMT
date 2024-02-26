using Plots # Must be installed with "add Plots"
using DelimitedFiles # No need to install

# Example providing path with strings (not always safe in windows)
base = "examples/data/example3/"
writedlm(base*"RPM_d0.15.csv", rpm_range) # rpm_range is not a vector but a UnitRange
writedlm(base*"T_d0.15.csv", T) # write the vector to file
writedlm(base*"Q_d0.15.csv", Q) # "*" used to combine strings

# Example providing path with joinpath
RPM_file = readdlm(joinpath(base,"RPM_d0.15.csv"))
T_file = readdlm(joinpath(base,"T_d0.15.csv"))
Q_file = readdlm(joinpath(base,"Q_d0.15.csv"))


# Plot results
p1 = plot(
    RPM_file, T_file, 
    label="Thrust", xlabel="RPM", ylabel="T [N]"
    )
p2 = plot(
    RPM_file, Q_file,
    label="Torque", xlabel="RPM", ylabel="Q [Nm]"
    )

# Plot with some plot attributes (see Plots.jl docs)
plot(
    p1,p2, 
    framestyle=:box,
    fg_legend=:false,
    layout=(2,1),
    size=(600,500)
    ) 

savefig(joinpath(examples_dir,"example3_loaded_results.svg"))