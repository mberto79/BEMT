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

thrust_element(cl, cd, c, U_corr, phi, ρ) = begin
    qA = 0.5*ρ*U_corr^2*c # (c*dr)
    ( cl*cos(phi) - cd*sin(phi) )*qA
end

thrust_momentum(vc, vi, r, ρ) = begin
    4*π*ρ*(vc + vi)*vi*r # (rdr)
end

element_momentum_balance(vc, vi, rpm, r, θ, cl, cd, chord) = begin
    Ω = (2π/60)*rpm
    U_r = Ω*r
    ϕ = atan((vc + vi)/U_r) 
    U_corr = sqrt( (vc + vi)^2 + (U_r)^2 )
    α = θ(r) - ϕ
    qA = 0.5*U_corr^2*chord # (ρ*dr) missing
    Te =( cl(α)*cos(ϕ) - cd(α)*sin(ϕ) )*qA
    Tm = 4*π*(vc + vi)*vi*r # (ρ*dr) missing
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
    func, funcValue; 
    guess_range::Tuple{T,T}=(), args=(), tol::T=1e-8, trials::I=10,show=false
    ) where {T<:Number,I<:Integer}
    
    f(x) = func(x, args...,) # redefine input function to include arguments

    g1, g2 = guess_range
    g = g1:g2/trials:g2
    fg = zeros(length(g))
    for i ∈ eachindex(fg)
        fg[i] = f(g[i]) - funcValue
    end

    idx = sortperm(fg)
    g1 = g[idx[1]]
    g2 = g[idx[2]]

    f1 = f(g1) - funcValue
    f2 = f(g2) - funcValue
    
    for i ∈ 1:20
    g2, g1, f1 = secant_method(g1, g2, f1, f2)
    f2 = f(g2) - funcValue
        if abs(f2) <= tol
            show ? println("Converged: ", i, " iterations.") : nothing
            return g2
        end
    end
    println("Warning: Convergence criterion not met!")
end


"""
xn, x2, f2 = secant_method(x1, x2, f1, f2)

Direct implementation of the secant method equation. Takes current guesses (`x1` amd `x2`) and their function evaluations (`f1` and `f2`), returning the next guess `xn` and the previous state (`x2`, `f2``)
"""
function secant_method(x1, x2, f1, f2)
    xn = x2 - f2*(x2 - x1)/(f2 - f1)
    return xn, x2, f2
end