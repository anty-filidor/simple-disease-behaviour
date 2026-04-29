"""
    jacob(x, z, R, p, c) -> Matrix

Return the 2×2 Jacobian of the system evaluated at state `(x, z)`.
"""
function jacob(x, z, R, p, c)
    J = [
        -p*(1 - 6x + 6x^2 + z - 2x*z) - c*(1-z),
        R*z*(1-z),
        -(1-x)*(p*x + c),
        R*x*(1-2z) - 1,
    ]
    J = reshape(J, 2, 2)
    return J
end

"""
    compute_eigenvalues_jacob(x, z, R, p, c) -> Tuple

Return eigenvalues and eigenvectors of the Jacobian at state `(x, z)`.
"""
function compute_eigenvalues_jacob(x, z, R, p, c)
    J = jacob(x, z, R, p, c)
    e = eigen(J)
    vals = e.values
    vectors = e.vectors
    return vals, vectors
end

"""
    calculate_values(R, p, c) -> Vector{Float64}

Return all x-coordinates of endemic equilibria in (0, 1) for the given parameters.
"""
function calculate_values(R, p, c)
    poly_coeff = [-c/p, -1, 2*R, -2*R]
    f(x) = poly_coeff[1] + poly_coeff[2] * x + poly_coeff[3] * x^2 + poly_coeff[4] * x^3
    return find_zeros(f, 0, 1)
end

"""
    find_dfe(R, p, c) -> Tuple

Return `(eigenvalues, eigenvectors, is_stable)` for the disease-free equilibrium (x=1, I=0).

The DFE is stable when R < 1 (eigenvalue R−1 < 0) and unstable otherwise.
"""
function find_dfe(R, p, c)
    vals, vectors = compute_eigenvalues_jacob(1.0, 0.0, R, p, c)
    is_stable = all(real.(vals) .<= 0)
    return vals, vectors, is_stable
end

"""
    find_equilibria(R, p, c) -> Vector

Return all endemic equilibria as `[x*, I*, is_stable]` for the given parameters.
"""
function find_equilibria(R, p, c)
    roots = calculate_values(R, p, c)
    equilibria = []
    for x in roots
        z = (R*x - 1) / (R*x)
        z > 0 || continue
        values, _ = compute_eigenvalues_jacob(x, z, R, p, c)
        real_values = real(values)
        if all(real_values .<= 0)
            push!(equilibria, (x, z, true))
        else
            push!(equilibria, (x, z, false))
        end
    end
    return equilibria
end
