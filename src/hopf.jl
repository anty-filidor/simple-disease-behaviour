"""
    hopf_data(x, z, R, d) -> Union{NamedTuple, Nothing}

For an endemic equilibrium `(x*, z*)`, return `(p_c, ω, det_J)` if a Hopf bifurcation exists
(i.e. `A ≠ 0`, `p_c > 0`, `det(J) > 0`), otherwise `nothing`.
"""
function hopf_data(x::Float64, z::Float64, R::Float64, d::Float64)
    A = -(1 - 6x + 6x^2 + z - 2x*z) - d*(1 - z)
    B = R*x*(1 - 2z) - 1
    abs(A) < 1e-12 && return nothing
    p_c = -B / A
    p_c > 0 || return nothing
    c_c = p_c * d
    J = jacob(x, z, R, p_c, c_c)
    det_J = det(J)
    det_J > 0 || return nothing
    ω = sqrt(det_J)
    return (p_c = p_c, ω = ω, det_J = det_J)
end

"""
    find_hopf_points(R, d) -> Vector

Return Hopf data for all non-trivial endemic equilibria at `(R, d)`.
Each entry contains `(x, z, p_c, ω, det_J)`.
"""
function find_hopf_points(R::Float64, d::Float64)
    xs = calculate_values(R, 1.0, d)
    results = []
    for x in xs
        z = (R*x - 1) / (R*x)
        z > 0 || continue
        hd = hopf_data(x, z, R, d)
        isnothing(hd) && continue
        push!(results, (x = x, z = z, p_c = hd.p_c, ω = hd.ω, det_J = hd.det_J))
    end
    return results
end
