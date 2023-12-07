export BGeometry, uniform_mesh
export load_aerofoil, sigma, thrust_coefficient, integrate
export linear_twist
export calculate_vi
export induced_angle, corrected_velocity
export thrust_element, thrust_momentum, trust_balance

struct BGeometry{I<:Integer,F<:AbstractFloat,V<:AbstractArray} 
    n_blades::I
    n_panels::I
    n_edges::I
    r::V
    radius::F
    dr::F
end

uniform_mesh(radius, n_blades, n_panels) = begin
    n_edges = n_panels+1
    dr = radius/n_panels
    r = [0.0:dr:radius;]
    return BGeometry(
        n_blades,
        n_panels,
        n_edges,
        r,
        radius,
        dr
    )    
end

load_aerofoil(file; startline=0) = begin
    data = readdlm(file, ',', Float64, skipstart=startline-1)
    alpha = deg2rad.(data[:,1])
    cl = Spline1D(alpha, data[:,2]) # generate cl fit
    cd = Spline1D(alpha, data[:,3]) # generate cd fit
    cl, cd
end

linear_twist(a_root, a_tip, radius) = begin
    k1 = deg2rad(a_root)
    k2 = deg2rad(a_tip)
    (r) -> k1 + (k2 - k1)/radius*r
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

induced_angle(vc, vi, rpm, r) = begin
    atan((vc + vi)/((2π/60)*rpm*r))
end

corrected_velocity(vc, vi, rpm, r) = begin
    sqrt( (vc + vi)^2 + ((2π/60)*rpm*r)^2 )
end

thrust_element(rotor, vi, vc, rpm, ρ, cl, cd, θ, chord) = begin
    (; r, n_blades, n_edges) = rotor
    dT = zeros(eltype(r), n_edges)
    Ω = (2π/60)*rpm
    for i ∈ eachindex(r)
        ri = r[i]
        vii = vi[i]
        U_r = Ω*ri
        ϕ = atan((vc + vii)/U_r) 
        U_corr = sqrt( (vc + vii)^2 + (U_r)^2 )
        α = θ(ri) - ϕ
        qA = 0.5*ρ*U_corr^2*chord(ri) # (dr) missing
        dT[i] =( cl(α)*cos(ϕ) - cd(α)*sin(ϕ) )*qA*n_blades
    end
    return dT
end

thrust_momentum(rotor, vi, vc, ρ) = begin
    r = rotor.r
    dT = zeros(eltype(r), rotor.n_edges)
    for i ∈ eachindex(r)
        dT[i] = 4*π*ρ*(vc + vi[i])*vi[i]*r[i]
    end
    return dT
end

trust_balance(vi, vc, Ω, nb, r, θ, chord, cl, cd) = begin
    U_r = Ω*r
    ϕ = atan((vc + vi)/U_r) 
    U_corr = sqrt( (vc + vi)^2 + (U_r)^2 )
    α = θ(r) - ϕ
    qA = 0.5*U_corr^2*chord(r) # (ρ*dr) missing
    Te =( cl(α)*cos(ϕ) - cd(α)*sin(ϕ) )*qA*nb
    Tm = 4*π*(vc + vi)*vi*r # (ρ*dr) missing
    Te - Tm
end

calculate_vi(rotor, vc, rpm, θ, chord, cl, cd) = begin
    (; r, radius, n_blades) = rotor
    Ω = (2π/60)*rpm
    v_tip = Ω*radius
    vi = similar(r)
    for i ∈ eachindex(r)
        args = (vc, Ω, n_blades, r[i], θ, chord, cl, cd)
        vi[i] = secant_solver(
            trust_balance, 0.0, guess_range=(0.0, v_tip/2), args=args)
    end
    return vi
end

integrate(fx, x) = begin
    n = length(fx)
    sum = zero(Float64)
    for i ∈ 1:(n-1)
        sum += 0.5*(fx[i] + fx[i+1])*(x[i+1] - x[i])
    end
    sum 
end