using Plots
using BEMT

file = "naca001264-il-200000-n5.csv"
cl, cd = load_aerofoil(file; startline=12)
a = [-0.22:0.001:0.22;]

plot(rad2deg.(a), cl.(a))
plot!(rad2deg.(a), cd.(a))

diameter = 25.40e-2
radius = diameter/2
chord = 0.15*radius
nb = 2
rpm = 3000
omega = rpm*(2π/60)
T = 1.82
v_tip = omega*radius
# J = V/(nD) # n in rev/sec
n = rpm*(1/60)
J = 0.6
vc = J*n*diameter
# vc = 100

ρ = 1.225
σ = sigma(chord, radius, nb)

induced_velocity(Vc, T, ρ, A) = -Vc/2 + sqrt(Vc/2^2 + T/(2*ρ*A))
area(r) = r^2*π
phi(r, vc, vi, rpm) = atan((vc + vi)/((2π/60)*rpm*r))
lambda(r, vc, vi, rpm) = phi(r, vc, vi, rpm)*r
slope(r) = 2π*r
pitch(r, at, ar) = deg2rad(at) - deg2rad(at - ar)*(r)

vi = induced_velocity(vc, T, ρ, area(radius))
r = [0.25:0.01:1.0;]
λ = lambda.(r, vc, vi, rpm)
ϕ = phi.(r,  vc, vi, rpm)
a = slope.(r)
θ = pitch.(r, 12, 10)
sig = solidity(r)

dCt = thrust_coefficient.(σ, a, θ, r, λ)
Ct = integrate(dCt, r)
T = integrate(dCt, r)*area(radius)*ρ*v_tip^2
CtUI = T/(ρ*n^2*(diameter)^4)

thrust.(vc, vi, rpm, r, θ, cl, cd, ρ) 
T_new = integrate(dT, r)*area(radius)*ρ*v_tip^2

# plot(r, dCt)
# scatter([J], [CtUI], label="$rpm")
scatter!([J], [CtUI], label=:none)
# scatter([vc], [T], label=:none)
# scatter!([vc], [T], label=:none)
# scatter([rpm], [T], label="$rpm")
# plot(r, rad2deg.(ϕ))
# plot(r, λ)
plot(r, dCt)
plot!(r, dCt_new)