# Parameter Space Analysis: `hopf_phase_Rd.png`

## The bug in the original figure

The original figure coloured each cell by `max tr(J)` across endemic equilibria.
This caused the entire small-d region to appear as deep stable blue, which was wrong.

**Why `tr(J) < 0` does not imply stability.** For a 2×2 system, stability requires
_both_ `tr(J) < 0` _and_ `det(J) > 0`. A saddle point has `det(J) < 0` — eigenvalues
of opposite sign — so one eigenvalue is always positive (unstable). But the trace
`tr = λ₁ + λ₂` can still be negative if the negative eigenvalue dominates in magnitude.

Concrete example at R = 5, p = 1.05, d = 0.05: the single endemic equilibrium is a
saddle with tr(J) ≈ −3.18 and eigenvalues λ ≈ +0.19, −3.37. The old figure showed
deep stable blue (tr ≈ −3). The corrected figure shows red (max Re(λ) = +0.19, unstable).

## The fix

The corrected diagnostic is: for each cell, compute `max Re(λ)` per endemic equilibrium,
then take the **minimum over all equilibria**.

- `min < 0` — at least one stable endemic equilibrium exists → **blue**
- `min > 0` — all equilibria are unstable → **red** (stable limit cycle)

This correctly catches saddles: even when tr < 0, max Re(λ) > 0 if det < 0.

## Parameter space structure

| Region        | Colour | What exists                        |
| ------------- | ------ | ---------------------------------- |
| Disease-free  | Grey   | No endemic equilibrium             |
| Saddle region | Red    | One endemic EE, a saddle (det < 0) |
| Stable region | Blue   | Stable focus + saddle              |
| Hopf-unstable | Red    | Unstable focus (tr > 0)            |

The two red regions together are the **stable limit cycle domain**: all endemic
equilibria are unstable, the system is bounded in [0,1]², so by Poincaré–Bendixson
a stable limit cycle must exist.

## The black contour

The contour at `min Re(λ_max) = 0` traces two boundaries:

- **Hopf bifurcation curve** (upper, bubble-shaped): equilibrium crosses from stable
  focus to unstable focus; limit cycle is born here.
- **Transcritical boundary** (lower, nearly horizontal): `d = 1/R − 2/R²`. Below
  this curve no stable focus exists yet; above it the stable focus appears alongside
  the saddle.

Together they enclose the full red (limit cycle) domain.
