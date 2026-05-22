"""
    hopf_data(x, z, R, d) -> Union{NamedTuple, Nothing}

For an endemic equilibrium `(x*, z*)`, return `(p_c, ω, det_J)` if a Hopf bifurcation exists
(i.e. `A ≠ 0`, `p_c > 0`, `det(J) > 0`), otherwise `nothing`.
"""
function hopf_data(x::Float64, z::Float64, R::Float64, d::Float64)
    # tr(J) = J[1,1] + J[2,2]
    #       = [-p(1-6x+6x²+z-2xz) - c(1-z)] + [Rx(1-2z) - 1]
    # Substitute c = p·d and collect by p:  tr(J) = p·A + B

    A = -(1 - 6x + 6x^2 + z - 2x*z) - d*(1 - z)  # coefficient of p in tr(J)
    B = R*x*(1 - 2z) - 1                            # p-independent part of tr(J) = J[2,2]

    abs(A) < 1e-12 && return nothing  # tr(J) independent of p → no unique Hopf crossing
    p_c = -B / A                      # tr(J) = 0  →  p·A = -B  →  p_c = -B/A
    p_c > 0 || return nothing         # negative imitation rate is unphysical

    c_c = p_c * d                     # recover c at the bifurcation point via c = p·d

    J = jacob(x, z, R, p_c, c_c)     # Jacobian evaluated at the Hopf point (where tr = 0)
    det_J = det(J)
    det_J > 0 || return nothing       # det(J) > 0 required: eigenvalues ±iω must be purely imaginary
    ω = sqrt(det_J)                   # at tr = 0 eigenvalues are ±i·√det(J), so ω = √det(J)
    return (p_c = p_c, ω = ω, det_J = det_J)
end

"""
    find_hopf_points(R, d) -> Vector

Return Hopf data for all non-trivial endemic equilibria at `(R, d)`.
Each entry contains `(x, z, p_c, ω, det_J)`.
"""
function find_hopf_points(R::Float64, d::Float64)
    # The equilibrium polynomial 2Rx³-2Rx²+x+d = 0 depends only on R and d = c/p,
    # not on p and c individually, so p=1, c=d is a valid proxy for any p at this d.
    xs = calculate_values(R, 1.0, d)
    results = []
    for x in xs
        z = (R*x - 1) / (R*x)    # I* from dI/dt = 0 with I > 0: R·x·(1-I) = 1 → I* = (Rx-1)/(Rx)
        z > 0 || continue          # skip roots where the implied I* ≤ 0 (no endemic equilibrium)
        hd = hopf_data(x, z, R, d)
        isnothing(hd) && continue  # no valid Hopf point at this equilibrium
        push!(results, (x = x, z = z, p_c = hd.p_c, ω = hd.ω, det_J = hd.det_J))
    end
    return results
end
