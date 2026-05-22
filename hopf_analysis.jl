using GLMakie
using OrdinaryDiffEq
using Roots
using LinearAlgebra

include("src/params.jl")
include("src/model.jl")
include("src/equilibria.jl")
include("src/hopf.jl")

d0 = c0 / p0   # d = c/p: ratio of spontaneous switching to imitation, used as the y-axis in Fig A
println("Default parameters: R₀=$R0, p₀=$p0, c₀=$c0, d₀=$(round(d0, digits=4))")

# ─── Figure A: (R, d) heatmap ────────────────────────────────────────────────
# Sweep (R, d) at fixed p = p₀ and colour each cell by max Re(λ) across all
# endemic equilibria, where λ are eigenvalues of J.
# max Re(λ) > 0 → unstable;  max Re(λ) < 0 → stable.
# Cells with no endemic equilibrium are left as NaN and rendered in grey.
#
# NOTE: tr(J) alone is NOT a valid stability diagnostic.  A saddle point has
# det(J) < 0 (eigenvalues of opposite sign), so one eigenvalue is positive even
# when tr(J) < 0.  Using max Re(λ) correctly identifies saddles as unstable.

NR, Nd = 400, 400
R_grid = range(2.0, 8.0, length = NR)   # reproduction number axis
d_grid = range(0.01, 0.30, length = Nd) # d = c/p axis

max_re_eig = fill(NaN32, NR, Nd)  # preallocate; NaN flags cells with no endemic equilibrium

println("Computing max Re(λ) heatmap on $(NR)×$(Nd) grid…")
for j in eachindex(d_grid)
    d = d_grid[j]
    c = p0 * d   # recover c from d = c/p at the fixed p₀
    for i in eachindex(R_grid)
        R = R_grid[i]
        xs = calculate_values(R, p0, c)   # x-coordinates of all endemic equilibria
        re_vals = Float64[]
        for x in xs
            z = (R*x - 1) / (R*x)   # corresponding I* from dI/dt = 0
            z > 0 || continue         # skip spurious roots with no positive prevalence
            J = jacob(x, z, R, p0, c)
            push!(re_vals, maximum(real.(eigvals(J))))   # max Re(λ): positive = unstable
        end
        isempty(re_vals) && continue                       # no endemic equilibrium at this (R, d)
        max_re_eig[i, j] = Float32(maximum(re_vals))      # take the most-unstable equilibrium per cell
    end
end

fig_a = Figure(size = (750, 550))
ax_a = Axis(
    fig_a[1, 1],
    xlabel = "R",
    ylabel = "d = C/p",
    title = "Max Re(λ) at endemic equilibria  (p = p₀ = $p0)",
)

hm = heatmap!(
    ax_a,
    R_grid,
    d_grid,
    max_re_eig,
    colormap = Reverse(:RdBu),
    colorrange = (-1, 1),
    nan_color = :lightgray,
)

# Re(λ) = 0 stability boundary — replace NaN with a large negative value so the
# contour algorithm never interpolates across the disease-free (grey) region.
stability_for_contour = replace(max_re_eig, NaN32 => -2.0f0)
contour!(
    ax_a,
    R_grid,
    d_grid,
    stability_for_contour,
    levels = [0.0],
    color = :black,
    linewidth = 2,
)

scatter!(
    ax_a,
    [R0],
    [d0],
    color = :limegreen,
    markersize = 14,
    marker = :star5,
    strokecolor = :black,
    strokewidth = 1,
    label = "(R₀, d₀)",
)

Colorbar(fig_a[1, 2], hm, label = "max Re(λ)")
axislegend(ax_a, position = :rb)

mkpath(figs_path)
save("$figs_path/hopf_phase_Rd.png", fig_a, px_per_unit = 1)
println("Saved $figs_path/hopf_phase_Rd.png")

# ─── Figure B: (R, p) Hopf curves ────────────────────────────────────────────
# Fix d = d₀ and sweep R.  For each R, find_hopf_points returns the critical p_c
# (and oscillation frequency ω) at every endemic equilibrium.  There can be up to
# two equilibria (EE₁, EE₂) — each traces its own branch in the (R, p) plane.

R_vals = collect(range(2.0, 8.0, length = 600))

# Separate accumulators for the first and second equilibrium branches.
branch1_R, branch1_pc, branch1_ω = Float64[], Float64[], Float64[]
branch2_R, branch2_pc, branch2_ω = Float64[], Float64[], Float64[]

println("Tracing Hopf curves in (R, p) plane at d = d₀…")
for R in R_vals
    pts = find_hopf_points(R, d0)   # returns 0, 1, or 2 Hopf points at this R
    if length(pts) >= 1             # first endemic equilibrium (EE₁)
        push!(branch1_R, R)
        push!(branch1_pc, pts[1].p_c)   # critical p at which EE₁ loses stability
        push!(branch1_ω, pts[1].ω)      # oscillation frequency at that point
    end
    if length(pts) >= 2             # second endemic equilibrium (EE₂), if it exists
        push!(branch2_R, R)
        push!(branch2_pc, pts[2].p_c)
        push!(branch2_ω, pts[2].ω)
    end
end

println(
    "EE₁ Hopf branch: $(length(branch1_R)) points, R ∈ [$(round(minimum(branch1_R), digits=2)), $(round(maximum(branch1_R), digits=2))]",
)
if !isempty(branch2_R)
    println(
        "EE₂ Hopf branch: $(length(branch2_R)) points, R ∈ [$(round(minimum(branch2_R), digits=2)), $(round(maximum(branch2_R), digits=2))]",
    )
end

fig_b = Figure(size = (750, 500))
ax_p = Axis(
    fig_b[1, 1],
    xlabel = "R",
    ylabel = "p",
    title = "Hopf bifurcation curves  (d = d₀ ≈ $(round(d0, digits=3)))",
)

!isempty(branch1_R) && lines!(
    ax_p,
    branch1_R,
    branch1_pc,
    color = :darkred,
    linewidth = 2,
    label = "EE₁  p_c(R)",
)
!isempty(branch2_R) && lines!(
    ax_p,
    branch2_R,
    branch2_pc,
    color = :darkorange,
    linewidth = 2,
    linestyle = :dash,
    label = "EE₂  p_c(R)",
)
hlines!(
    ax_p,
    [p0],
    color = :black,
    linewidth = 1.5,
    linestyle = :dot,
    label = "p = p₀ = $p0",
)
axislegend(ax_p, position = :rt)
xlims!(ax_p, 2.5, 4.5)
ylims!(ax_p, 0.0, p0 * 6)

ax_ω = Axis(
    fig_b[1, 1],
    yaxisposition = :right,
    ylabel = "ω",
    yticklabelcolor = :dodgerblue,
    ytickcolor = :dodgerblue,
    rightspinecolor = :dodgerblue,
    backgroundcolor = :transparent,
)
hidespines!(ax_ω, :l, :t, :b)
hidexdecorations!(ax_ω)
linkxaxes!(ax_p, ax_ω)

!isempty(branch1_R) && lines!(
    ax_ω,
    branch1_R,
    branch1_ω,
    color = :dodgerblue,
    linewidth = 1.5,
    linestyle = :dot,
    label = "EE₁  ω(R)",
)
!isempty(branch2_R) && lines!(
    ax_ω,
    branch2_R,
    branch2_ω,
    color = :steelblue,
    linewidth = 1.5,
    linestyle = :dashdot,
    label = "EE₂  ω(R)",
)

save("$figs_path/hopf_phase_Rp.png", fig_b, px_per_unit = 1)
println("Saved $figs_path/hopf_phase_Rp.png")
