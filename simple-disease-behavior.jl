
import Pkg
Pkg.add("DynamicalSystems")
Pkg.add("OrdinaryDiffEq")
Pkg.add("GLMakie")
Pkg.add("Roots")

using DynamicalSystems
using OrdinaryDiffEq
using GLMakie
using Base.Threads
using Random
using Roots
using LinearAlgebra

T = 100
Ttr = 0
sampling_time = 1
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

# check what eigen vectors are in equilibria

R0 = 3.0
p0 = 1.05
c0 = 0.19

N = 2

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

Equilibria = find_equilibria(R0, p0, c0)

function cond(x, I)
    if x >= 0 && x <= 1 && I >= 0 && I <= 1 && I > 0 
        return true
    else
        return false
    end
end

function dynamics_rule_si(x0, param, t)
    R, p, c = param
    eps = 0.0
    f = 15
    R_sin = R + eps * sin(2 * π * t / f)
    x, I = x0
    dxdt = -p*(1+I)*x*(1-x)^2 + p*(1-I)*(1-x)*x^2 + c*(1-I)*(1-x)
    dIdt = R_sin*x*I*(1-I) - I
    r = SVector{2}([dxdt; dIdt])
    return r
end

function my_trajectory(x0, I, r, p, c)
    if cond(x0, I)
        X0 = [x0 I]
        pars = [r p c]
        dynamics = ContinuousDynamicalSystem(dynamics_rule_si, X0, pars; diffeq)
        SOL, t = trajectory(dynamics, T; Ttr=Ttr, Δt=sampling_time)
        t0 = Int((Ttr + 1) / sampling_time)
        x = SOL[t0:end, 1]
        I = SOL[t0:end, 2]
        t = t[t0:end]
        println(size(I))
        println(size(t))
        return [x, I]
    end
end

function my_trajectory_no_condition(x0, I, r, p, c)
    X0 = [x0 I]
    pars = [r p c]
    dynamics = ContinuousDynamicalSystem(dynamics_rule_si, X0, pars; diffeq)
    SOL, t = trajectory(dynamics, T; Δt=sampling_time, Ttr=Ttr)
    x = SOL[:, 1]
    I = SOL[:, 2]
    println(size(I))
    return [x, I]
end

# P = vec(collect(Iterators.product(LinRange(0.0, 0.15, N), LinRange(0.1, 0.3, N), LinRange(0.25, 0.35, N), LinRange(0.0, 0.75, N))))

# P = vec(collect(Iterators.product(LinRange(0.5, 1.0, N), LinRange(0.5, 1.0, N))))
P = vec(collect(Iterators.product([0.4], [0.25])))
P = vec(collect(Iterators.product([0.4], [0.2])))
# P = [[0.6, 0.3]]


# unit_v = [Equilibria[1][1] - Equilibria[2][1],
#             Equilibria[1][2] - Equilibria[2][2],
#             Equilibria[1][3] - Equilibria[2][3],
#             Equilibria[1][4] - Equilibria[2][4]]
# println("Unit vector: $unit_v")
# Δ = sqrt(sum(unit_v .* unit_v)) / 150
# println("Δ: $Δ")
# unit_v = unit_v ./ sqrt(sum(unit_v .* unit_v))
# EEb_eq = Equilibria[2][1:4]
# P = [EEb_eq .+ unit_v.*Δ*i for i in 148:150]

# P = rand(Point3f, 10000)

# Compute all posible trajectories
data = Array{Any}(undef, length(P))
@threads for ind in eachindex(P)
    x = Float64(P[ind][1])
    I = Float64(P[ind][2])
    res = my_trajectory(x, I, R0, p0, c0)
    if !isnothing(res)
        data[ind] = res
    end
end


f2 = Figure(size = (600, 450))

ax = Axis(f2[1, 1],
    title = "Experiment R = $R0, p = $p0, c=$c0",
    ylabel = "I",
    xlabel = "x",
    # title = "",
    # xlabel = "",
    # ylabel = "",
    # zlabel = "",
    # aspect=(1.0, 1.0, 1.0),
)


for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:end-1]).^2 .+ (I[2:end] - I[1:end-1]).^2)
    x = x[1:end-1]
    I = I[1:end-1]
    lines!(ax,
            x,
            I,
            colormap=cgrad(Makie.to_colormap(:Reds)[4:9]),
            alpha=0.3,
            color=cs,
            overdraw = true,
            # label="DFE1 trajectories"
            )
end

println(Equilibria)
if length(Equilibria) > 0
    scatter!(ax,
        Equilibria[1][1],
        Equilibria[1][2],
        color=:red,
        marker = Equilibria[1][3] == true ? :circle : :star4,
        markersize=20.0,
        alpha=1.0,
        strokecolor=:white,
        strokewidth=1,
        label="EE1")
end


if length(Equilibria) > 1
    scatter!(ax,
        Equilibria[2][1],
        Equilibria[2][2],
        color=:blue,
        marker = Equilibria[1][3] == true ? :circle : :star4,
        markersize=20.0,
        alpha=1.0,
        strokecolor=:white,
        strokewidth=1,
        label="EE1")
end

out_path = "simple-disease-behavior-plots"
mkpath(out_path)
save(f"./{out_path}/f1.png", f2, px_per_unit=1)
display(f2)

ft = Figure(size = (600, 450))

ax = Axis(ft[1, 1],
    title = "Experiment R = $R0, p = $p0, c=$c0",
    ylabel = "I",
    xlabel = "t",
    # title = "",
    # xlabel = "",
    # ylabel = "",
    # zlabel = "",
    # aspect=(1.0, 1.0, 1.0),
)


for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:end-1]).^2 .+ (I[2:end] - I[1:end-1]).^2)
    x = x[1:end-1]
    I = I[1:end-1]
    lines!(ax,
            I,
            colormap=cgrad(Makie.to_colormap(:Reds)[4:9]),
            alpha=0.3,
            color=cs,
            overdraw = true,
            # label="DFE1 trajectories"
            )
end

save(f"./{out_path}/ft.png", ft, px_per_unit=1)
display(f2)
# ft

ft = Figure(size = (600, 450))

ax = Axis(ft[1, 1],
    title = "Experiment R = $R0, p = $p0, c=$c0",
    ylabel = "I",
    xlabel = "t",
    # title = "",
    # xlabel = "",
    # ylabel = "",
    # zlabel = "",
    # aspect=(1.0, 1.0, 1.0),
)


for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:end-1]).^2 .+ (I[2:end] - I[1:end-1]).^2)
    x = x[1:end-1]
    I = I[1:end-1]
    lines!(ax,
            I,
            colormap=cgrad(Makie.to_colormap(:Reds)[4:9]),
            alpha=0.3,
            color=cs,
            overdraw = true,
            # label="DFE1 trajectories"
            )
end
# display(ft)
