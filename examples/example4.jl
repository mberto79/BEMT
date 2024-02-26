using Plots # Must be installed with "add Plots"
using DelimitedFiles # No need to install
using Printf # provides macro to format strings 

d_string = @sprintf "%.2f" 25e-2 # useful macro to convert numbers to formatted strings

# Example providing path with strings (not always safe in windows)
base = "examples/data/example4/"
writedlm(base*"RPM_d$(d_string).csv", rpm_range) # use "$" to interpolate "d_string"
writedlm(base*"T_d$(d_string).csv", T) 
writedlm(base*"Q_d$(d_string).csv", Q) 

# Example providing path with joinpath
RPM_file = readdlm(joinpath(base,"RPM_d$(d_string).csv"))
T_file = readdlm(joinpath(base,"T_d$(d_string).csv"))
Q_file = readdlm(joinpath(base,"Q_d$(d_string).csv"))
