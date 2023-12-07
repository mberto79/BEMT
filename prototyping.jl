using Plots
using BEMT

file = "naca001264-il-200000-n5.csv"
cl, cd = load_aerofoil(file; startline=12)
a = [-0.22:0.001:0.22;]

# plot(rad2deg.(a), cl.(a))
# plot(rad2deg.(a), cd.(a))

diameter = 25.40e-2
radius = diameter/2

nb = 2
rpm = 1000
omega = rpm*(2π/60)
v_tip = omega*radius
# J = V/(nD) # n in rev/sec
n = rpm*(1/60)
J = 0.0
vc = J*n*diameter
# vc = 0
ρ = 1.225


n = 75

rotor = uniform_mesh(radius, nb, n)
θ = linear_twist(20, 7.5, radius)
chord(r) = 0.15*radius

vi = similar(rotor.r_edges)
@time for i ∈ eachindex(rotor.r_edges)
    args = (vc, rpm, nb, rotor.r_edges[i], θ, cl, cd, chord)
    vi[i] = secant_solver(
        trust_balance, 0.0, guess_range=(0.0, v_tip/2), args=args)
end


dTm = thrust_momentum(rotor, vi, vc, ρ)
sum(dTm)
thrust = integrate(dTm, rotor.r_edges)

plot(rotor.r_edges, dTm)
plot!(rotor.r_edges, vi)