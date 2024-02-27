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
diameter = 15e-2
radius = diameter/2
nb = 3 # number of blades

# Define mesh for BET
n = 100
rotor = uniform_mesh(radius, nb, n) # returns a "meshed" rotor object

# Define geometry (these are functions)
theta = linear_function(deg2rad(20), deg2rad(15), radius) # 20 to 7.5 degree twist
# chord = constant_function(0.035) # constant chord of 0.075
# chord = linear_function(0.05, 0.02, rotor.radius) # linear taper from 0.05 to 0.02 m
# chord = nonlinear_function([0.0, 0.25, 0.5, 1.0].*radius, [0.05, 0.04, 0.04, 0.025]) 

# Operating conditions
vc = 0
rho = 1.225
# rmp = 500 # Original setting for a single RPM run
rpm_range = 100:100:5000 # We want to loop over the RPM variable

# Pre-allocate output (single values are now vectors)

T = zeros(length(rpm_range))
Q = zeros(length(rpm_range))
P = zeros(length(rpm_range))
dT_all = zeros(elements(rotor),length(rpm_range))
last_iter = length(rpm_range)

for (i, rpm) in enumerate(rpm_range)
    println("Evaluating rpm: ", rpm)

    omega = rpm*(2π/60)
    v_tip = omega*radius

    # Solve BEMT equations to determine induced velocity
    vi, converged = calculate_vi(rotor, vc, rpm, theta, chord, cl, cd, warnings=false)
    
    if !converged
        last_iter = i-1
        println("Last converged iteration ", last_iter)
        break 
    end

    # Calculate aerodynamic performance
    dT, dQ, dP = element_performance(rotor, vi, vc, rpm, rho, cl, cd, theta, chord)
    dT_all[:,i] .= dT # Store thrust distribution along span (notice "." for performance)

    # Integrate element results over the rotor blades
    T[i] = integrate1(dT, rotor.r) # Rotor thrust prediction (by BEM)
    Q[i] = integrate1(dQ, rotor.r) # Rotor torque prediction
    P[i] = integrate1(dP, rotor.r) # Rotor power prediction

end

# Plot results
p1 = plot(
    rpm_range[1:last_iter], T[1:last_iter], 
    label="Thrust", xlabel="RPM", ylabel="T [N]"
    )
p2 = plot(
    rpm_range[1:last_iter], Q[1:last_iter], 
    label="Torque", xlabel="RPM", ylabel="Q [Nm]"
    )
p3 = plot(
    rpm_range[1:last_iter], P[1:last_iter], 
    label="Power", xlabel="RPM", ylabel="Q [W]"
    )

p4 = plot()
for (i, rpm) in enumerate(rpm_range)
    plot!(
        p4, rotor.r, dT_all[:,i], 
        xlabel="Radius [m]", ylabel="Thurst [N/m]", legend=false)
end

rpm_range[1:10:end]

plot(p1,p2,p3,p4) 

# savefig(joinpath(examples_dir,"example2_results.svg"))