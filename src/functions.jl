export sigma, thrust_coefficient 

sigma(chord, radius, nb) = begin
    nb*chord/(π*radius)
end

thrust_coefficient(r, σ, a, θ, λ) = begin
    n = length(r)
    dT = zeros(n)
    Ct = zero(Float64)
    for i ∈ 1:(n-1)
        r1 = r[i]
        r2 = r[i+1]
        h = r2 - r1
        f1 = (σ/2)*a[i]  *(θ[i]  *r1^2 - λ[i]  *r[i]  )
        f2 = (σ/2)*a[i+1]*(θ[i+1]*r2^2 - λ[i+1]*r[i+1])
        Ct += 0.5*(f1 + f2)*h
        dT[i]   = (σ/2)*a[i]  *(θ[i]  *r1^2 - λ[i]  *r[i]  )
        dT[i+1] = (σ/2)*a[i+1]*(θ[i+1]*r2^2 - λ[i+1]*r[i+1])
    end
    Ct, dT 
end

