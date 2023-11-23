export sigma, thrust_coefficient, integrate

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

