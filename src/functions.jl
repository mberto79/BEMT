export BGeometry, radius, elements, points, uniform_mesh
export load_xfoil, load_xflr5, sigma, integrate
export nonlinear_function, linear_function, constant_function
export calculate_vi
export element_performance, thrust_momentum, trust_balance

struct BGeometry{I<:Integer,F<:AbstractFloat,V<:AbstractArray} 
    n_blades::I # rotor blades
    n_panels::I # number of panels/elements in mesh
    n_edges::I # number of edges in mesh
    r::V # vector of coordinate along radius
    radius::F # radius of rotor
    dr::F # space between mesh points
end

# BGeometry access functions
radius(rotor::BGeometry) = rotor.r
elements(rotor::BGeometry) = rotor.n_panels
points(rotor::BGeometry) = rotor.n_edges

uniform_mesh(radius, n_blades, n_panels) = begin
    n_edges = n_panels+1
    dr = radius/(n_panels-1)
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

load_xfoil(file; startline=12) = begin
    data = readdlm(file, ',', Float64, skipstart=startline-1)
    cl, cd = fit_polar(data)
    return cl, cd
end

load_xflr5(file; startline=12) = begin
    lines = readlines(file)
    nlines = length(lines)
    data_lines = nlines-startline-1
    out = zeros(data_lines, 3)
    line_count = 0
    for i ∈ startline:nlines
        line = lines[i]
        sline = split(line)
        line_count += 1
        data_count = 0 
        for item ∈ sline 
            data = tryparse(Float64, item)
            if data !== nothing && data_count <=2
                data_count +=1
                out[line_count, data_count] = data
            end
        end
    end
    cl, cd = fit_polar(out)
    return cl, cd
end

fit_polar(data) = begin 
    alpha = deg2rad.(data[:,1])
    cl = Spline1D(alpha, data[:,2]) # generate cl fit
    cd = Spline1D(alpha, data[:,3]) # generate cd fit
    cl, cd
end

nonlinear_function(x_vals, y_vals) = begin
    Spline1D(x_vals, y_vals)
end

linear_function(val_root, val_tip, radius) = begin
    (r) -> val_root + (val_tip - val_root)/radius*r
end

constant_function(k) = begin
    (r) -> k
end

sigma(chord, radius, nb) = begin
    nb*chord/(π*radius)
end

element_performance(rotor, vi, vc, rpm, ρ, cl, cd, θ, chord) = begin
    (; r, n_blades, n_panels, n_edges, dr) = rotor
    dT = similar(r)
    dQ = similar(r)
    dP = similar(r)
    Ω = (2π/60)*rpm
    for i ∈ eachindex(r)
        ri = r[i]
        vii = vi[i]
        U_r = Ω*ri
        ϕ = atan((vc + vii)/U_r) 
        U_corr = sqrt( (vc + vii)^2 + (U_r)^2 )
        α = θ(ri) - ϕ
        qN = 0.5*ρ*U_corr^2*chord(ri)*n_blades # (dr) missing
        sinϕ = sin(ϕ); cosϕ = cos(ϕ)
        clα = cl(α); cdα = cd(α)
        dT[i] = ( clα*cosϕ - cdα*sinϕ )*qN
        dQ[i] = ( clα*sinϕ + cdα*cosϕ )*qN*ri
        dP[i] = ( clα*sinϕ + cdα*cosϕ )*qN*ri*Ω
    end
    return dT, dQ, dP
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

calculate_vi(
    rotor, vc, rpm, θ, chord, cl, cd; show_convergence=false, warnings=true
    ) = begin
    (; r, radius, n_blades) = rotor
    Ω = (2π/60)*rpm
    v_tip = Ω*radius
    vi = similar(r); converged = false
    for i ∈ eachindex(r)
        args = (vc, Ω, n_blades, r[i], θ, chord, cl, cd)
        vi[i], converged = secant_solver(
            trust_balance, 0.0, guess_range=(0.0, v_tip*2), args=args, show_convergence=show_convergence, warnings=warnings
            )
        if !converged
            vi[i] = 0.0
            break
        end
    end
    return vi, converged
end

integrate(fx, x) = begin
    n = length(fx)
    sum = zero(Float64)
    for i ∈ 1:(n-1)
        sum += (fx[i] + fx[i+1])*(x[i+1] - x[i])
    end
    sum 
end