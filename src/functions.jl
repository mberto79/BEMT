export load_aerofoil, sigma, thrust_coefficient, integrate, thrust

load_aerofoil(file; startline=0) = begin
    data = readdlm(file, ',', Float64, skipstart=startline-1)
    alpha = deg2rad.(data[:,1])
    cl = Spline1D(alpha, data[:,2]) # generate cl fit
    cd = Spline1D(alpha, data[:,3]) # generate cd fit
    cl, cd
end

sigma(chord, radius, nb) = begin
    nb*chord/(π*radius)
end

# Induced angle
phi(r, vi, vc, rpm) = atan((vc + vi)/((2π/60)*rpm*r))

# True velocity 
U_inf(r, vi, vc, rpm) = sqrt( (vc + vi)^2 + ((2π/60)*rpm*r)^2 )

# Lift and drag forces 
aerofoil_forces(ρ, U_inf, cl, cd, c, dr) = begin
    lift = 0.5*ρ*U_inf^2*c*cl*dr
    drag = 0.5*ρ*U_inf^2*c*cd*dr
    lift, drag
end

thrust_coefficient(σ, a, θ, r, λ) = begin
    (σ/2)*a*(θ*r^2 - λ*r)
end

get_thrust(cl, cd, σ, ϕ, r) = begin
    (σ/2)*(cl*cos(ϕ) - cd*sin(ϕ))*r^2
end

thrust(vc, vi, rpm, r, theta, cl, cd, ρ) = begin
    phi = atan((vc + vi)/((2π/60)*rpm*r))
    alpha = theta - phi
    U_corr = sqrt( (vc + vi)^2 + ((2π/60)*rpm*r)^2 )
    qA = 0.5*ρ*U_corr^2*c # (c*dr)
    lift = qA*cl(alpha)
    drag = qA*cd(alpha)
    lift*cos(phi) - drag*sin(phi)
end

integrate(fx, x) = begin
    n = length(fx)
    sum = zero(Float64)
    for i ∈ 1:(n-1)
        sum += 0.5*(fx[i] + fx[i+1])*(x[i+1] - x[i])
    end
    sum 
end