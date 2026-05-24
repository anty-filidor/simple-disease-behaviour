# Parameter Space Analysis: `hopf_phase_Rd.png`

## What prompted the investigation

The original figure (using `max tr(J)`) showed two visually distinct blue regions
separated by a sharp border at roughly d ≈ 0.12. This was initially dismissed as
a colourmap saturation artefact (tr hitting the −1 clip). That was wrong.

## What the border actually is

The border is the curve **d = 1/R − 2/R²**. This is where the lower endemic
equilibrium emerges with I\* → 0, i.e. it crosses the disease-free boundary into
biological validity. Below this curve the second equilibrium does not yet exist
(its implied I\* < 0); above it a stable focus appears alongside the existing saddle.

## What lived on each side of the border in the old figure

| Old colour | d relative to border | Endemic equilibria            | Eigenvalues                                        |
| ---------- | -------------------- | ----------------------------- | -------------------------------------------------- |
| Deep navy  | below                | 1 — a **saddle**              | λ ≈ +0.19, −3.37 (tr ≈ −3.18, det ≈ −0.65)         |
| Light blue | above, outside Hopf  | 2 — **stable focus** + saddle | focus: λ ≈ −0.13 ± 0.32i; saddle: λ ≈ +0.20, −3.29 |
| Red        | inside Hopf bubble   | 1 — unstable focus            | tr > 0                                             |

The deep navy region was appearing stable because tr(J) ≈ −3. It is not stable.

## The bug: `tr(J) < 0` does not imply stability

Stability of a 2×2 equilibrium requires **both** tr(J) < 0 **and** det(J) > 0.
A saddle has det(J) < 0 — eigenvalues of opposite sign — so one is always positive
(the equilibrium is unstable). But tr = λ₁ + λ₂ can be negative if the negative
eigenvalue dominates, which is exactly what happened here (λ ≈ +0.19, −3.37 → tr ≈ −3).

The original `max tr(J)` heatmap therefore could not distinguish a saddle from a
stable node, and showed the entire saddle region (small d) as deep stable blue.

## The fix: `min over EE of max Re(λ)`

For each grid cell: compute `max Re(λ)` per endemic equilibrium, then take the
minimum over all equilibria.

- `min < 0` — at least one stable EE exists → **blue**
- `min > 0` — all EEs are unstable → **red** (stable limit cycle by Poincaré–Bendixson)

This correctly marks saddles as unstable regardless of their trace.

## Parameter space structure (corrected figure)

| Region                         | Colour | Dynamics                             |
| ------------------------------ | ------ | ------------------------------------ |
| Disease-free                   | Grey   | I → 0                                |
| Saddle region (d < 1/R − 2/R²) | Red    | One saddle EE; no stable fixed point |
| Stable region                  | Blue   | Stable focus + saddle coexist        |
| Hopf-unstable (inside bubble)  | Red    | Unstable focus; stable limit cycle   |

The full red area is the **stable limit cycle domain**.

## The black contour

The contour at `min Re(λ_max) = 0` now traces two boundaries that together
enclose the entire limit cycle domain:

- **Hopf bifurcation curve** (upper, bubble-shaped): limit cycle is born here.
- **Transcritical boundary** (lower, d = 1/R − 2/R²): stable focus is born here.
