using Plots
using BEMT

file = "naca001264-il-200000-n5.csv"
cl, cd = load_aerofoil(file; startline=12)
a = [-0.22:0.001:0.22;]

# plot(rad2deg.(a), cl.(a))
# plot(rad2deg.(a), cd.(a))

diameter = 25.40e-2
radius = diameter/2
chord = 0.15*radius
nb = 2
rpm = 3000
omega = rpm*(2π/60)
v_tip = omega*radius
# J = V/(nD) # n in rev/sec
n = rpm*(1/60)
J = 0.6
vc = J*n*diameter
# vc = 0
ρ = 1.225


r = radius.*[0.25:0.01:1.0;]
θ = linear_twist(15, 2)

vi = similar(r)

blade_vi(vi, vc, rpm, r, θ, cl, cd, chord) = begin
    element_momentum_balance(vc, vi, rpm, r, θ, cl, cd, chord)
end

@time for i ∈ eachindex(r)
    args = (vc, rpm, r[i], θ, cl, cd, chord)
    vi[i] = secant_solver(
        blade_vi, 0.0, guess_range=(0.0, v_tip/2), args=args)
end


alphai = induced_angle.(vc, vi, rpm, r) 
U_corr = corrected_velocity.(vc, vi, rpm, r)
alpha = θ.(r) .- alphai
dT = thrust_element.(cl.(alpha), cd.(alpha), chord, U_corr, alphai, ρ)
dTm = thrust_momentum.(vc, vi, r, ρ)
thrust = integrate(dT, r)

plot(r, dT)
plot!(r, dTm)
plot!(r, vi)

scatter!([J], [thrust], label="J = $J")

a = (1,2)

g1, g2 = a