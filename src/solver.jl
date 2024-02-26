export secant_solver

"""
function secant_solver(
    func, funcValue; 
    guess_range::Tuple{T,T}=(), 
    args=(), 
    tol::T=1e-8, 
    trials::I=10,
    show=false
    ) where {T<:Number,I<:Integer}

Find the value x for which a function `func` will return the target value `funcValue`, given lower and upper guesses, `g1` and `g2` to an absolute tolerance `tol` (set to 1e-8 as default)
"""
function secant_solver(
    func, funcValue; 
    guess_range::Tuple{T,T}=(), args=(), tol::T=1e-8, trials::I=10,
    show_convergence=false, warnings=true
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
    
    for i ∈ 1:200
    g2, g1, f1 = secant_method(g1, g2, f1, f2)
    f2 = f(g2) - funcValue
        if abs(f2) <= tol
            show_convergence ? println("Converged: ", i, " iterations.") : nothing
            return g2, true
        end
    end

    if warnings
        print("\nWarning: Convergence criterion not met!\n")
        print("Check 1: Input configuration is not physical (reduce/increase RPM or θ)\n")
        print("Check 2: Angle operating beyond range used for cl and cd (lower θ or Vc)\n")
        print("Check 3: Mesh does not have sufficient points\n")
        print("Check 4: Solver tolerance too high\n\n")
    end

    return 0.0, false
end


"""
xn, x2, f2 = secant_method(x1, x2, f1, f2)

Direct implementation of the secant method equation. Takes current guesses (`x1` amd `x2`) and their function evaluations (`f1` and `f2`), returning the next guess `xn` and the previous state (`x2`, `f2``)
"""
function secant_method(x1, x2, f1, f2)
    xn = x2 - f2*(x2 - x1)/(f2 - f1)
    return xn, x2, f2
end