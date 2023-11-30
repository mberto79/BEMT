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
T = 3.5
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
θ = pitch.(r, 5, 2)
# sig = solidity(r)

dCt = thrust_coefficient.(σ, a, θ, r, λ)
Ct = integrate(dCt, r)
T = integrate(dCt, r)*area(radius)*ρ*v_tip^2
CtUI = T/(ρ*n^2*(diameter)^4)

U_corr = zeros(length(r))
alpha = zeros(length(r))
alphai = zeros(length(r))

alphai .= induced_angle.(vc, vi, rpm, r) 
U_corr .= corrected_velocity.(vc, vi, rpm, r)
alpha .= θ .- alphai
dT = thrust_element.(cl.(alpha), cd.(alpha), chord, U_corr, alphai, ρ)
dTm = thrust_momentum.(vc, vi, r, ρ)
T_new = integrate(dT, r)

plot(r, dT)
plot!(r, dTm)