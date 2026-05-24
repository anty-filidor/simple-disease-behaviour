using GLMakie
using OrdinaryDiffEq
using Roots
using LinearAlgebra

include("src/params.jl")
include("src/model.jl")
include("src/equilibria.jl")
include("src/hopf.jl")
include("src/hopf_plots.jl")

d0 = c0 / p0   # d = c/p: ratio of spontaneous switching to imitation, used as the y-axis in Fig A
println("Default parameters: R₀=$R0, p₀=$p0, c₀=$c0, d₀=$(round(d0, digits=4))")

mkpath(figs_path)

plot_heatmap_Rd(R0, p0, d0, figs_path)   # Fig A: (R, d) stability heatmap
plot_hopf_Rp(p0, d0, figs_path)          # Fig B: Hopf critical curves in (R, p) plane
plot_limit_cycle(R0, p0, c0, figs_path)  # Fig C: phase portrait of the stable limit cycle
