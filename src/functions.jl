export load_aerofoil, sigma, thrust_coefficient, integrate, get_thrust

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

integrate(fx, x) = begin
    n = length(fx)
    sum = zero(Float64)
    for i ∈ 1:(n-1)
        sum += 0.5*(fx[i] + fx[i+1])*(x[i+1] - x[i])
    end
    sum 
end

thrust_coefficient(σ, a, θ, r, λ) = begin
    (σ/2)*a*(θ*r^2 - λ*r)
end

get_thrust(cl, cd, σ, ϕ, r) = begin
    (σ/2)*(cl*cos(ϕ) - cd*sin(ϕ))*r^2
end

thrust(cl, cd, r) = begin
    ϕ = atan((vc + vi)/((2π/60)*rpm*r))
    u_inf = sqrt( (rpm*y*2π/60)^2 + (Vc + Vi)^2 )
    (cl*cos(ϕ) - cd*sin(ϕ))*(0.5*ρ*u_inf^2*area)
end
