using DynamicalSystems
using OrdinaryDiffEq
using GLMakie
using Base.Threads
using Random
using Roots
using LinearAlgebra

include("src/params.jl")
include("src/model.jl")
include("src/equilibria.jl")
include("src/trajectories.jl")

# --- Equilibria ---

Equilibria = find_equilibria(R0, p0, c0)
dfe_vals, _, dfe_stable = find_dfe(R0, p0, c0)
println("DFE stable: $dfe_stable, eigenvalues: $dfe_vals")

# --- Initial conditions ---

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

# --- Simulate trajectories ---

data = Array{Any}(undef, length(P))
@threads for ind in eachindex(P)
    x = Float64(P[ind][1])
    I = Float64(P[ind][2])
    res = my_trajectory(x, I, R0, p0, c0)
    if !isnothing(res)
        data[ind] = res
    end
end

# --- Phase portrait (x vs I) ---

f2 = Figure(size = (600, 450))

ax =
    Axis(f2[1, 1], title = "Experiment R = $R0, p = $p0, c=$c0", ylabel = "I", xlabel = "x")

for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:(end-1)]) .^ 2 .+ (I[2:end] - I[1:(end-1)]) .^ 2)
    x = x[1:(end-1)]
    I = I[1:(end-1)]
    lines!(
        ax,
        x,
        I,
        colormap = cgrad(Makie.to_colormap(:Reds)[4:9]),
        alpha = 0.3,
        color = cs,
        overdraw = true,
    )
end

println(Equilibria)
if length(Equilibria) > 0
    scatter!(
        ax,
        Equilibria[1][1],
        Equilibria[1][2],
        color = :red,
        marker = Equilibria[1][3] == true ? :circle : :star4,
        markersize = 20.0,
        alpha = 1.0,
        strokecolor = :white,
        strokewidth = 1,
        label = "EE1",
    )
end

if length(Equilibria) > 1
    scatter!(
        ax,
        Equilibria[2][1],
        Equilibria[2][2],
        color = :blue,
        marker = Equilibria[1][3] == true ? :circle : :star4,
        markersize = 20.0,
        alpha = 1.0,
        strokecolor = :white,
        strokewidth = 1,
        label = "EE1",
    )
end

mkpath(figs_path)
save("$figs_path/f1.png", f2, px_per_unit = 1)
display(f2)

# --- Time series (I vs t) ---

ft = Figure(size = (600, 450))

ax =
    Axis(ft[1, 1], title = "Experiment R = $R0, p = $p0, c=$c0", ylabel = "I", xlabel = "t")

for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:(end-1)]) .^ 2 .+ (I[2:end] - I[1:(end-1)]) .^ 2)
    x = x[1:(end-1)]
    I = I[1:(end-1)]
    lines!(
        ax,
        I,
        colormap = cgrad(Makie.to_colormap(:Reds)[4:9]),
        alpha = 0.3,
        color = cs,
        overdraw = true,
    )
end

save("$figs_path/ft.png", ft, px_per_unit = 1)
display(f2)
# ft

ft = Figure(size = (600, 450))

ax =
    Axis(ft[1, 1], title = "Experiment R = $R0, p = $p0, c=$c0", ylabel = "I", xlabel = "t")

for ind in filter(i -> isassigned(data, i), 1:length(data))
    x = data[ind][1]
    I = data[ind][2]
    cs = sqrt.((x[2:end] - x[1:(end-1)]) .^ 2 .+ (I[2:end] - I[1:(end-1)]) .^ 2)
    x = x[1:(end-1)]
    I = I[1:(end-1)]
    lines!(
        ax,
        I,
        colormap = cgrad(Makie.to_colormap(:Reds)[4:9]),
        alpha = 0.3,
        color = cs,
        overdraw = true,
    )
end
# display(ft)
