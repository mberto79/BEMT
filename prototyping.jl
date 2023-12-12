using Plots
using BEMT

file_xfoil = "naca001264-il-200000-n5.csv"
cl, cd = load_xfoil(file_xfoil)

file_xflr = "NACA 0012_T1_Re0.030_M0.00_N4.0.txt"
cl, cd = load_xflr5(file_xflr)

a = [-0.22:0.001:0.22;]

plot(rad2deg.(a), cl.(a))
plot(rad2deg.(a), cd.(a))

diameter = 25.40e-2
radius = diameter/2
ρ = 1.225
nb = 2
rpm = 2500
omega = rpm*(2π/60)
v_tip = omega*radius

# J = V/(nD) # n in rev/sec
# n = rpm*(1/60)
# J = 0.0
# vc = J*n*diameter
vc = 0

n = 75
rotor = uniform_mesh(radius, nb, n)
θ = linear_twist(20, 7.5, radius)
chord(r) = 0.15*radius

vi = calculate_vi(rotor, vc, rpm, θ, chord, cl, cd)

dT, dQ, dP = element_performance(rotor, vi, vc, rpm, ρ, cl, cd, θ, chord)
dTm = thrust_momentum(rotor, vi, vc, ρ)
thrust = integrate(dTm, rotor.r)
thrust = integrate(dT, rotor.r)
torque = integrate(dQ, rotor.r)
power = integrate(dP, rotor.r)

plot(rotor.r, dT)
plot(rotor.r, dQ)
plot!(rotor.r, dP)
plot!(rotor.r, vi)


scatter!(
    [rpm], [torque], 
    label="RPM = $rpm",
    xlabel="RPM", ylabel="Torque [N]"
    )

