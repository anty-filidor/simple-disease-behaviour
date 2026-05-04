# Hopf Bifurcation Analysis

## Context

Section 5 of the paper is currently empty. The goal is to characterise, as explicitly as possible,
the parameter domain in which the non-trivial endemic equilibria undergo a Hopf bifurcation — i.e.
the boundary between a stable endemic equilibrium and a stable limit cycle.

The key analytical fact: at an endemic equilibrium `(x*, z*)` the Jacobian trace splits as

    tr(J) = p · A(x*, z*, d)  +  B(x*, z*, R)

where `d = C/p`, `A = -(1 - 6x + 6x² + z - 2xz) - d(1-z)`, and `B = Rx(1-2z) - 1`.
Setting `tr(J) = 0` gives the explicit critical behavioural strength:

    p_c(R, d)  =  -B / A       (provided A ≠ 0 and det(J) > 0)

At the Hopf point the eigenvalues are purely imaginary `±iω` with `ω = √(det J)|_{p=p_c}`,
giving the limit-cycle frequency directly.

Because the system can have two non-trivial endemic equilibria (EE₁, EE₂), there are
potentially two Hopf curves in each parameter plane. Both must be tracked.

---

## Steps

### 0. Checkout branch and Bump version

Checkout to a new branch: `feature--hpf-bifurcation-analysis`

Increment the patch version in `Project.toml` from `0.2.1` to `0.3.0` (minor bump — new
analysis output).

---

### 1. Add `src/hopf.jl`

Create a new source file with the Hopf utility functions. It depends only on `equilibria.jl`
(already included before it in any script).

Key functions to implement (signatures only — bodies are for the dev-assistant):

```julia
"""
    hopf_data(x, z, R, d) -> Union{NamedTuple, Nothing}

For an endemic equilibrium `(x*, z*)`, return `(p_c, ω, det_J)` if a Hopf bifurcation exists
(i.e. `A ≠ 0`, `p_c > 0`, `det(J) > 0`), otherwise `nothing`.
"""
function hopf_data(x::Float64, z::Float64, R::Float64, d::Float64)

"""
    find_hopf_points(R, d) -> Vector{NamedTuple}

Return Hopf data for all non-trivial endemic equilibria at `(R, d)`.
Each entry contains `(x, z, p_c, ω, det_J)`.
"""
function find_hopf_points(R::Float64, d::Float64)
```

`find_hopf_points` must call `calculate_values(R, 1.0, d)` (using dummy `p=1`, `c=d`, since
the cubic only depends on `d = c/p`) and filter for `z > 0`.

---

### 2. Extend `test/runtests.jl`

Add a `"Hopf bifurcation"` testset that verifies correctness of the new functions:

- **Algebraic identity**: at `(x*, z*, R, p_c, p_c·d)`, check `tr(J) ≈ 0` and `det(J) > 0`.
- **Sign flip**: at `p_c ± ε`, check that `tr(J)` changes sign as expected.
- **Frequency consistency**: check that `Im(λ) ≈ ω` and `Re(λ) ≈ 0` for the eigenvalues of
  `jacob(x*, z*, R, p_c, p_c·d)`.
- **No Hopf without endemic equilibrium**: for `R = 1.5` (below invasion threshold), check
  that `find_hopf_points` returns an empty vector.
- **Default parameters cross-check**: at `(R₀=3.0, d₀=c₀/p₀≈0.181)`, check that `p_c` is
  finite, positive, and that `p₀=1.05` sits in the unstable (limit-cycle) region, i.e.
  `tr(J) > 0` at the endemic equilibrium with default params.

Run with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

### 3. Write `hopf_analysis.jl`

A standalone script (analogous to `main.jl`) that includes the same source files and produces
two figures saved to `figs/`.

**Figure A — `figs/hopf_phase_Rd.png`**: `(R, d)` parameter plane at fixed `p = p₀`.

- Grid: `R ∈ [2, 8]`, `d ∈ [0.01, 0.30]`, ~400 × 400 points.
- For each grid point, evaluate `tr(J)` at all non-trivial endemic equilibria (using `jacob`
  with actual `p = p₀`, `c = p₀·d`). Store the maximum `tr` value.
- Heatmap coloured by that value (`RdBu` diverging palette, clamped to `[-1, 1]`).
- Black contour at `tr = 0` — the Hopf curve. Two contour lines if two EE branches are
  present in that region.
- Grey / `NaN` where no non-trivial endemic equilibrium exists.
- Mark the default operating point `(R₀, d₀)` as a scatter point.

**Figure B — `figs/hopf_phase_Rp.png`**: `(R, p)` parameter plane at fixed `d = d₀ = c₀/p₀`.

- For each `R ∈ [2, 8]`, call `find_hopf_points(R, d₀)` to get `p_c` for each EE branch.
- Plot `p_c(R)` as explicit curves (one per branch, labelled EE₁ / EE₂).
- Add a horizontal dashed line at `p = p₀` to show where the default parameters sit.
- On a twin y-axis (or inset), overlay the Hopf frequency `ω(R)` along each branch.

Both figures should have axis labels and titles that match the paper's notation (`R`, `d = C/p`,
`p`).

---

### 4. Verify visually and numerically

After running `julia --project=. hopf_analysis.jl`:

- **Consistency with existing trajectories**: the default point `(R₀, p₀, c₀)` should lie in
  the limit-cycle region on both figures (it is known from `main.jl` that the trajectory
  oscillates at default params).
- **Boundary sanity**: the Hopf curve should not appear in the grey "no endemic equilibrium"
  region.
- **Two-equilibrium regime**: for `R > 4d₀ + 2 ≈ 2.72`, confirm that two EE branches are
  visible and each has its own `p_c` curve on Figure B.
- **Trajectory spot-check**: pick one point clearly inside the predicted limit-cycle region and
  one clearly outside; run `main.jl` with those params and confirm the trajectories match the
  prediction (oscillation vs. convergence to a fixed point).
