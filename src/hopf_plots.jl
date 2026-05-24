"""
    plot_heatmap_Rd(R0, p0, d0, figs_path)

Figure A: (R, d) parameter plane coloured by `min over EE of max Re(λ)`.
Blue = stable fixed point exists; red = all EEs unstable (limit cycle domain);
grey = no endemic equilibrium.  Black contour is the stability boundary.
"""
function plot_heatmap_Rd(R0, p0, d0, figs_path)
    NR, Nd = 400, 400
    R_grid = range(2.0, 8.0, length = NR)
    d_grid = range(0.01, 0.30, length = Nd)

    stability_map = fill(NaN32, NR, Nd)   # NaN → grey (no endemic EE)

    println("Computing stability heatmap on $(NR)×$(Nd) grid…")
    for j in eachindex(d_grid)
        d = d_grid[j]
        c = p0 * d                             # recover c from d = c/p at fixed p₀
        for i in eachindex(R_grid)
            R = R_grid[i]
            xs = calculate_values(R, p0, c)
            re_vals = Float64[]
            for x in xs
                z = (R*x - 1) / (R*x)
                z > 0 || continue
                J = jacob(x, z, R, p0, c)
                push!(re_vals, maximum(real.(eigvals(J))))  # max Re(λ) per EE; >0 ↔ unstable
            end
            isempty(re_vals) && continue
            stability_map[i, j] = Float32(minimum(re_vals))  # most-stable EE
        end
    end

    fig = Figure(size = (750, 550))
    ax = Axis(
        fig[1, 1],
        xlabel = "R",
        ylabel = "d = c/p",
        title = "Stability of most-stable EE  (p = p₀ = $p0)",
    )

    hm = heatmap!(
        ax,
        R_grid,
        d_grid,
        stability_map,
        colormap = Reverse(:RdBu),
        colorrange = (-1, 1),
        nan_color = :lightgray,
    )

    # NaN → −2 so the contour algorithm never interpolates across the grey region
    stability_for_contour = replace(stability_map, NaN32 => -2.0f0)
    contour!(
        ax,
        R_grid,
        d_grid,
        stability_for_contour,
        levels = [0.0],
        color = :black,
        linewidth = 2,
    )

    scatter!(
        ax,
        [R0],
        [d0],
        color = :limegreen,
        markersize = 14,
        marker = :star5,
        strokecolor = :black,
        strokewidth = 1,
        label = "(R₀, d₀)",
    )

    Colorbar(fig[1, 2], hm, label = "min Re(λ_max)")
    axislegend(ax, position = :rb)

    save("$figs_path/hopf_phase_Rd.png", fig, px_per_unit = 1)
    println("Saved $figs_path/hopf_phase_Rd.png")
end

"""
    plot_hopf_Rp(p0, d0, figs_path)

Figure B: Hopf critical curve p_c(R) and oscillation frequency ω(R) at fixed d = d₀.
Separate branches are traced for EE₁ and EE₂ (when it exists).
"""
function plot_hopf_Rp(p0, d0, figs_path)
    R_vals = collect(range(2.0, 8.0, length = 600))

    branch1_R, branch1_pc, branch1_ω = Float64[], Float64[], Float64[]
    branch2_R, branch2_pc, branch2_ω = Float64[], Float64[], Float64[]

    println("Tracing Hopf curves in (R, p) plane at d = d₀…")
    for R in R_vals
        pts = find_hopf_points(R, d0)
        if length(pts) >= 1
            push!(branch1_R, R);
            push!(branch1_pc, pts[1].p_c);
            push!(branch1_ω, pts[1].ω)
        end
        if length(pts) >= 2
            push!(branch2_R, R);
            push!(branch2_pc, pts[2].p_c);
            push!(branch2_ω, pts[2].ω)
        end
    end

    println(
        "EE₁ Hopf branch: $(length(branch1_R)) points, " *
        "R ∈ [$(round(minimum(branch1_R), digits=2)), $(round(maximum(branch1_R), digits=2))]",
    )
    if !isempty(branch2_R)
        println(
            "EE₂ Hopf branch: $(length(branch2_R)) points, " *
            "R ∈ [$(round(minimum(branch2_R), digits=2)), $(round(maximum(branch2_R), digits=2))]",
        )
    end

    fig = Figure(size = (750, 500))
    ax_p = Axis(
        fig[1, 1],
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
        fig[1, 1],
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

    save("$figs_path/hopf_phase_Rp.png", fig, px_per_unit = 1)
    println("Saved $figs_path/hopf_phase_Rp.png")
end

"""
    plot_limit_cycle(R0, p0, c0, figs_path)

Figure C: phase portrait of the stable limit cycle.
Two trajectories — one started inside (near the unstable EE), one outside — converge
to the same orbit after discarding the transient.  Equilibria are marked with ✕ (unstable)
or ● (stable).
"""
function plot_limit_cycle(R0, p0, c0, figs_path)
    function ode!(du, u, p, _)   # in-place form; mirrors dynamics_rule_si without SVector
        R, im, c = p
        x, I = u
        du[1] = -im*(1+I)*x*(1-x)^2 + im*(1-I)*(1-x)*x^2 + c*(1-I)*(1-x)
        du[2] = R*x*I*(1-I) - I
    end

    function trajectory(x0, I0, T)
        prob = ODEProblem(ode!, [x0, I0], (0.0, T), [R0, p0, c0])
        sol = solve(prob, Vern9(), abstol = 1e-9, reltol = 1e-9, saveat = 0.05)
        return sol[1, :], sol[2, :]
    end

    eq = find_equilibria(R0, p0, c0)
    println("Endemic equilibria at (R₀=$R0, p₀=$p0, c₀=$c0):")
    for (i, e) in enumerate(eq)
        println(
            "  EE$i: x*=$(round(e[1], digits=4)), I*=$(round(e[2], digits=4)), stable=$(e[3])",
        )
    end

    # Outer start: perturb a cycle point radially outward by 20%.
    # A corner equilibrium at (1, 1-1/R) coexists with the limit cycle, so only
    # starting points inside the cycle's basin reach the attractor. Computing the
    # start from the settled orbit guarantees it lies just outside the cycle and
    # within that basin, regardless of parameter values.
    xs_settled, Is_settled = let
        prob = ODEProblem(
            ode!,
            [eq[1][1] + 0.001, eq[1][2] + 0.001],
            (0.0, 600.0),
            [R0, p0, c0],
        )
        sol = solve(prob, Vern9(), abstol = 1e-9, reltol = 1e-9, saveat = 0.05)
        idx = findfirst(≥(500.0), sol.t)
        sol[1, idx:end], sol[2, idx:end]
    end
    cx, cI = sum(xs_settled) / length(xs_settled), sum(Is_settled) / length(Is_settled)
    i_far = argmax(xs_settled)   # rightmost point on settled cycle
    x0_out = cx + 1.2 * (xs_settled[i_far] - cx)
    I0_out = cI + 1.2 * (Is_settled[i_far] - cI)

    x_in, I_in = trajectory(eq[1][1] + 0.02, eq[1][2] + 0.02, 100.0)   # near EE₁, inside
    x_out, I_out = trajectory(x0_out, I0_out, 100.0)                      # just outside the cycle

    fig = Figure(size = (600, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "x  (risky behaviour fraction)",
        ylabel = "I  (infected fraction)",
        title = "Stable limit cycle  (R₀=$R0, p₀=$p0, c₀=$c0)",
    )

    lines!(ax, x_in, I_in, color = :steelblue, linewidth = 2, label = "from inside")
    lines!(ax, x_out, I_out, color = :tomato, linewidth = 2, label = "from outside")

    for (i, e) in enumerate(eq)
        mk = e[3] ? :circle : :xcross   # ● stable, ✕ unstable
        scatter!(
            ax,
            [e[1]],
            [e[2]],
            marker = mk,
            markersize = 14,
            color = :black,
            label = "EE$i $(e[3] ? "(stable)" : "(unstable)")",
        )
    end

    axislegend(ax, position = :rt)
    save("$figs_path/hopf_limit_cycle.png", fig, px_per_unit = 1)
    println("Saved $figs_path/hopf_limit_cycle.png")
end
