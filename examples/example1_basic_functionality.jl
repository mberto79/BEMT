using Plots
using BEMT

# Define base directory where example datafiles are stored
examples_dir = pkgdir(BEMT, "examples")
data_dir = joinpath(examples_dir, "data")

# Define path for xfoil data file and loading
file_xfoil = joinpath(data_dir, "example1", "NACA001264_xfoil.csv")
cl, cd = load_xfoil(file_xfoil) # note cl and cd are functions

# Define path for fxlr5 data file and loading
file_xflr = joinpath(data_dir, "example1", "NACA0012_xflr5.txt")
cl, cd = load_xflr5(file_xflr) # note cl and cd are functions

# Define angle of attack range in degrees 
alpha = [-45:0.1:45;]

# Convert to radians
alpha_rad = deg2rad.(alpha) # apply deg2rad function element-wise (using ".")

# Plot alpha vs cl and cl 
plot(alpha, cl.(alpha_rad), label="Lift coefficient")

# Plot cd plot on top of previous plot - note "!" 
plot!(alpha, cd.(alpha_rad), label="Drag coefficient") 

# Add labels (documentation and options at Plots.jl)
plot!(xlabel="alpha", ylabel="Aerodynamic coefficients")


# Define propeller 
diameter = 25.40e-2
radius = diameter/2
nb = 2 # number of blades

# Operating conditions
vc = 0
rho = 1.225
rpm = 500
omega = rpm*(2π/60)
v_tip = omega*radius

# Define mesh for BET
n = 75
rotor = uniform_mesh(radius, nb, n)

# Define geometry (these are functions)
theta = linear_function(deg2rad(20), deg2rad(7.5), radius) # 20 to 7.5 degree twist
chord = constant_function(0.075) # constant chord of 0.075
# chord = linear_function(0.05, 0.02) # linear taper from 0.05 to 0.02 m

# Solve BEMT equations to determine induced velocity
vi, converged = calculate_vi(rotor, vc, rpm, theta, chord, cl, cd)

# Calculate aerodynamic performance
dT, dQ, dP = element_performance(rotor, vi, vc, rpm, rho, cl, cd, theta, chord)
dTm = thrust_momentum(rotor, vi, vc, rho) # used here as a check

# Integrate element results over the rotor blades
T = integrate(dT, rotor.r) # Rotor thrust prediction (by BEM)
Tm = integrate(dTm, rotor.r) # should be similar to the value of T
Q = integrate(dQ, rotor.r) # Rotor torque prediction
P = integrate(dP, rotor.r) # Rotor power prediction

# Plot results
p4 = plot(rotor.r, vi, label=:false, xlabel="Radius [m]", ylabel="Induced velocity")
p1 = plot(rotor.r, dT, label=:false, xlabel="Radius [m]", ylabel="Thrust / span")
p2 = plot(rotor.r, dQ, label=:false, xlabel="Radius [m]", ylabel="Torque / span")
p3 = plot(rotor.r, dP, label=:false, xlabel="Radius [m]", ylabel="Power / span")

plot(p1,p2,p3,p4, plot_title="Rotor Performance")
# savefig(joinpath(examples_dir,"example1_results.svg"))