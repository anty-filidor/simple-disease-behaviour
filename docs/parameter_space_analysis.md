# Parameter Space Analysis: `hopf_phase_Rd.png`

This document explains the structure of the (R, d) parameter plane, the stability
diagnostic used to colour it, and a bug in an earlier version of the figure.

## Parameter space structure

The heatmap sweeps the basic reproduction number R (x-axis) and the ratio
d = c/p (y-axis) at fixed p = p₀. Four qualitatively distinct regions exist:

| Region        | Colour | Condition                      | Dynamics                   |
| ------------- | ------ | ------------------------------ | -------------------------- |
| Disease-free  | Grey   | No endemic equilibrium         | I → 0                      |
| Saddle region | Red    | One endemic EE, a saddle       | No stable fixed point      |
| Stable region | Blue   | Two EEs: stable focus + saddle | Stable endemic fixed point |
| Hopf-unstable | Red    | One EE, unstable focus         | Stable limit cycle         |

The red regions (saddle + Hopf-unstable) together form the **stable limit cycle
domain**: when all endemic equilibria are unstable and the system is bounded in
[0,1]², the Poincaré–Bendixson theorem guarantees a stable limit cycle.

## The black contour

The contour at `min Re(λ_max) = 0` traces two distinct boundaries:

- **Hopf bifurcation curve** (upper, curved portion around the bubble): the endemic
  equilibrium transitions from stable focus to unstable focus. A limit cycle is born here.
- **Transcritical boundary** (lower, nearly horizontal curve): given by
  `d = 1/R − 2/R²` (requires R > 2). Below this curve the second endemic equilibrium
  does not exist (its implied I\* < 0); above it a stable focus appears alongside
  the saddle.

Together these two curves form the closed boundary of the entire limit cycle domain.

## The diagnostic: `min over EE of max Re(λ)`

For each grid cell, the code:

1. Finds all endemic equilibria (roots of the cubic equilibrium polynomial).
2. For each equilibrium, computes the maximum real part of the Jacobian eigenvalues: `max Re(λ)`.
3. Takes the **minimum** over all equilibria.

The result is negative if and only if at least one stable endemic equilibrium exists.

### Why `tr(J)` alone is insufficient

An earlier version of the figure used `max tr(J)` (maximum trace across equilibria).
This produced a misleading deep navy region at small d.

The error: `tr(J) < 0` does **not** imply stability. For a 2×2 Jacobian:

- Stability requires **both** `tr(J) < 0` **and** `det(J) > 0`.
- A **saddle point** has `det(J) < 0`, meaning one eigenvalue is positive and one is
  negative. The trace `tr = λ₁ + λ₂` can still be negative if the negative eigenvalue
  dominates in magnitude — even though the equilibrium is unstable.

**Numerical example** at R = 5, p = 1.05:

| d    | Equilibrium      | tr(J) | det(J) | max Re(λ) | Type              |
| ---- | ---------------- | ----- | ------ | --------- | ----------------- |
| 0.05 | upper (x ≈ 0.88) | −3.18 | −0.65  | **+0.19** | Saddle (unstable) |
| 0.15 | lower (x ≈ 0.22) | −0.25 | +0.12  | **−0.13** | Stable focus      |
| 0.15 | upper (x ≈ 0.86) | −3.08 | −0.67  | **+0.20** | Saddle (unstable) |

At d = 0.05 the old figure showed tr ≈ −3 → deep navy (appearing stable). The
corrected figure shows max Re(λ) = +0.19 → red (correctly unstable).

## Saddle region dynamics

Below the transcritical boundary (small d), only a single endemic equilibrium exists
and it is a saddle. Together with an unstable disease-free equilibrium (R > 1), the
phase portrait has no stable fixed point. By Poincaré–Bendixson the flow must
approach a limit cycle.

This region was entirely invisible in the original `tr(J)` figure — it looked
identical to the genuinely stable region.
