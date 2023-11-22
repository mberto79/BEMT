using Plots
using BEMT

chord = 0.1
radius = 2
nb = 3
σ = sigma(chord, radius, nb)

phi(r, vc, vi, rpm, R) = atan((vc + vi)/(0.10472*rpm*r*R))
lambda(r, vc, vi, rpm, R) = phi(r, vc, vi, rpm, R)*r
slope(r) = 2π
pitch(r, at, ar) = deg2rad(at) - deg2rad(at - ar)*(r)

r = [0.0:0.01:1.0;]
λ = lambda.(r, 0.0, 0.2, 500, 2)
ϕ = phi.(r,  0.0, 0.2, 500, 2)
a = slope.(r)
θ = pitch.(r, 10, 0)


Ct, dT = thrust_coefficient(r, σ, a, θ, λ)

plot(r, rad2deg.(ϕ))
plot(r, rad2deg.(θ))
plot(r, λ)
plot(r, dT)
Ct