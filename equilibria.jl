function jacob(x, z)
    R = R0
    p = p0
    c = c0
    J = [
        2*p*(1-x)*(2*x-(1+z)/2) - 2*p*(x*(x-(1+z)/2) + c*(1-z) / (2*p)),
        -p*(1+x)*(x+c/p),
        R*z*(1-z),
        R*x*(1-2*z) - 1
    ]
    J = reshape(J, 2, 2)
    return J
end

function compute_eigenvalues_jacob(x, z)
    J = jacob(x, z)
    e = eigen(J)
    vals = e.values
    vectors = e.vectors
    return vals, vectors
end

function calculate_values(R, p, c)
    # Define the polynomial coefficients
    poly_coeff = [
        -c/p,
        -1,
        2*R,
        -2*R
    ]
    f(x) = poly_coeff[1] + poly_coeff[2] * x + poly_coeff[3] * x^2 + poly_coeff[4] * x^3
    # Solve the polynomial
    roots = find_zeros(f, 0, 1)
    # Filter real roots
    real_roots = real(roots[abs.(imag.(roots)) .< 1e-17])
    return real_roots
end

function is_valid_solution(x, R, p, c)

    # Check for NaN or Inf values
    if any(isnan.([x])) || any(isinf.([x]))
        return false
    end

    # Check positivity
    if x < 0 || x > 1
        return false
    end

    return true
end

function find_equilibria(R, p, c)
    roots = calculate_values(R, p, c)
    equilibria = []
    for x in roots
        if is_valid_solution(x, R, p, c)
            z = (R*x-1)/(R*x)
            values, vectors = compute_eigenvalues_jacob(x, z)
            println("x: $x, z: $z, Eigenvalues: $values, Eigenvector: $vectors")
            real_values = real(values)
            if all(real_values .<= 0)
                push!(equilibria, [x, z, true])
            else
                push!(equilibria, [x, z, false])
            end
        end
    end
    return equilibria
end
