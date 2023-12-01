export load_aerofoil, sigma, thrust_coefficient, integrate
export linear_twist
export induced_angle, corrected_velocity
export thrust_element, thrust_momentum, element_momentum_balance
export secant_solver

load_aerofoil(file; startline=0) = begin
    data = readdlm(file, ',', Float64, skipstart=startline-1)
    alpha = deg2rad.(data[:,1])
    cl = Spline1D(alpha, data[:,2]) # generate cl fit
    cd = Spline1D(alpha, data[:,3]) # generate cd fit
    cl, cd
end

linear_twist(r, a_root, a_tip) = begin
    deg2rad(a_root) - deg2rad(a_root - a_tip)*(r)
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

thrust_element(cl, cd, c, U_corr, phi, ρ) = begin
    qA = 0.5*ρ*U_corr^2*c # (c*dr)
    ( cl*cos(phi) - cd*sin(phi) )*qA
end

thrust_momentum(vc, vi, r, ρ) = begin
    4*π*ρ*(vc + vi)*vi*r # (rdr)
end

element_momentum_balance(vi, vc, rpm, r, θ, cl, cd, chord, ρ) = begin
    ϕ = induced_angle(vc, vi, rpm, r) 
    U_corr = corrected_velocity(vc, vi, rpm, r)
    α = θ - ϕ
    Te = thrust_element(cl(α), cd(α), chord, U_corr, ϕ, ρ)
    Tm = thrust_momentum(vc, vi, r, ρ)
    Te - Tm
end

integrate(fx, x) = begin
    n = length(fx)
    sum = zero(Float64)
    for i ∈ 1:(n-1)
        sum += 0.5*(fx[i] + fx[i+1])*(x[i+1] - x[i])
    end
    sum 
end

"""
secant_solver(func, funcValue, g1, g2; tol=1e-8)

Find the value x for which a function `func` will return the target value `funcValue`, given lower and upper guesses, `g1` and `g2` to an absolute tolerance `tol` (set to 1e-8 as default)
"""
function secant_solver(
    func, funcValue, g1, g2; tol=1e-8, args=(), show=false)
    
    g = [g1:g2/10:g2;]
    fg = zeros(length(g))
    for i ∈ eachindex(fg)
        fg[i] = func(g[i], args...,) - funcValue
    end
    # fg = func.(g, args...,) .- funcValue
    I = sortperm(fg)
    g1 = g[I[1]]
    g2 = g[I[2]]

    f1 = func(g1, args...,) - funcValue
    f2 = func(g2, args...,) - funcValue
    
    for i ∈ 1:20
    g2, g1, f1 = secant_method(g1, g2, f1, f2)
    f2 = func(g2, args...,) - funcValue
        if abs(f2) <= tol
            show ? println("Converged: ", i, " iterations.") : nothing
            return g2
        end
    end
    println("Convergence criterion not met!")
end


"""
xn, x2, f2 = secant_method(x1, x2, f1, f2)

Direct implementation of the secant method equation. Takes current guesses (`x1` amd `x2`) and their function evaluations (`f1` and `f2`), returning the next guess `xn` and the previous state (`x2`, `f2``)
"""
function secant_method(x1, x2, f1, f2)
    xn = x2 - f2*(x2 - x1)/(f2 - f1)
    return xn, x2, f2
end