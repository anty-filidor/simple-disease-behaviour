# Fix plan

Addresses issues identified in the assessment of the coupled SIS + behaviour model.

---

## Commits

1. **Mathematical correctness** — Jacobian fix, parameterise `jacob`, remove dead seasonal forcing, remove no-op root filter
2. **DFE analysis** — add disease-free equilibrium finder and stability check
3. **Plot correctness** — fix equilibrium marker bug, fix time-series x-axis
4. **Dead code removal + docs** — delete unused function, debug prints, commented-out code, obsolete script; update README variable description

---

## Step 1 — Write tests first

File: `test/runtests.jl`

Setup: include `src/params.jl`, `src/model.jl`, `src/equilibria.jl`; use `Test` (stdlib) and `LinearAlgebra`. Add `Test` under `[extras]` and `[targets]` in `Project.toml`.

### `@testset "ODE RHS"`

| Test                               | How                                                          | Expected on current code |
| ---------------------------------- | ------------------------------------------------------------ | ------------------------ |
| `dIdt ≈ 0` at computed equilibrium | Evaluate `dynamics_rule_si` at `(x*, I*)` for default params | PASS                     |
| `dxdt ≈ 0` at computed equilibrium | Same call, check first component                             | PASS                     |
| `dIdt < 0` when `I > I*`           | Spot-check with `I = I* + 0.1`                               | PASS                     |

### `@testset "Boundary validator"`

| Test                         | How               | Expected |
| ---------------------------- | ----------------- | -------- |
| `cond(0.5, 0.3)` is `true`   | Direct call       | PASS     |
| `cond(-0.1, 0.3)` is `false` | x < 0             | PASS     |
| `cond(0.5, 0.0)` is `false`  | I = 0 is excluded | PASS     |
| `cond(0.5, 1.1)` is `false`  | I > 1             | PASS     |

### `@testset "Equilibrium finder"`

| Test                                           | How                             | Expected |
| ---------------------------------------------- | ------------------------------- | -------- |
| At least one root in (0, 1) for default params | `calculate_values(R0, p0, c0)`  | PASS     |
| Each root satisfies `I* = (R·x* − 1)/(R·x*)`   | Cross-check closed-form formula | PASS     |
| `dxdt ≈ 0` at each `(x*, I*)`                  | Plug into ODE                   | PASS     |
| `dIdt ≈ 0` at each `(x*, I*)`                  | Plug into ODE                   | PASS     |

### `@testset "Jacobian correctness"` — regression anchor for commit 1

Evaluate at fixed point `(x=0.4, I=0.2)`. Compare each element against its symbolic partial derivative.

| Element             | Symbolic formula               | Expected on current code              |
| ------------------- | ------------------------------ | ------------------------------------- |
| J[1,1] — ∂(dxdt)/∂x | derived from `dxdt` expression | PASS                                  |
| J[1,2] — ∂(dxdt)/∂I | `-(1-x)*(p*x + c)`             | **FAIL** — currently holds ∂(dIdt)/∂x |
| J[2,1] — ∂(dIdt)/∂x | `R*I*(1-I)`                    | **FAIL** — currently holds ∂(dxdt)/∂I |
| J[2,2] — ∂(dIdt)/∂I | `R*x*(1-2I) - 1`               | PASS                                  |

This testset must go red → green after commit 1.

### `@testset "Stability classification"`

| Test                                                | How                                       | Expected |
| --------------------------------------------------- | ----------------------------------------- | -------- |
| All eigenvalue real parts ≤ 0 at stable equilibrium | `compute_eigenvalues_jacob` at `(x*, I*)` | PASS     |
| Eigenvalue sign flips when R pushed below threshold | Use `R = 1.0`                             | PASS     |

Tolerances: `atol=1e-10` throughout (solver tolerance is `1e-9`).

---

## Commit 1 — Mathematical correctness

Files: `src/model.jl`, `src/equilibria.jl`

### A1 + A2 — Fix Jacobian and parameterise `jacob`

`src/equilibria.jl`

- Add `R, p, c` parameters to `jacob`; remove global reads on lines 3–5
- Reorder the four literal array elements so that after `reshape(J, 2, 2)` (column-major) the entries satisfy:
  - J[1,1] = ∂(dxdt)/∂x
  - J[1,2] = ∂(dxdt)/∂I = `-(1-x)*(p*x + c)`
  - J[2,1] = ∂(dIdt)/∂x = `R*I*(1-I)`
  - J[2,2] = ∂(dIdt)/∂I = `R*x*(1-2*I) - 1`
- Update `compute_eigenvalues_jacob` to pass `R, p, c` through to `jacob`
- Update `find_equilibria` to pass `R, p, c` through to `compute_eigenvalues_jacob`

Risk: high — verify each element symbolically before and after. The Jacobian testset must go green.

### A5 — Remove dead seasonal forcing

`src/model.jl`

- Delete `eps = 0.0` and `R_sin = R + eps * sin(2 * π * t / f)` and `f = 15`
- Replace `R_sin` with `R` (the parameter) in the `dIdt` line
- Keep the `t` argument in the function signature (required by DynamicalSystems.jl)

### A6 — Remove no-op polynomial root filter

`src/equilibria.jl:35`

- Replace `real(roots[abs.(imag.(roots)) .< 1e-17])` with `roots`
- `find_zeros` already returns `Float64`; the filter is always true

---

## Commit 2 — DFE analysis

File: `src/equilibria.jl`, `main.jl`

### A3 — Add `find_dfe`

`src/equilibria.jl`

- Add function `find_dfe(R, p, c)`:
  - DFE is at I = 0; x* satisfies dxdt = 0 with I = 0, i.e. solve `p*(1-x)_x^2 + c_(1-x) = 0` → x = 1 (and x = 0 for c = 0, but ignore)
  - Evaluate `jacob(x=1, I=0, R, p, c)` and return eigenvalues
  - Return whether DFE is stable

`main.jl`

- Call `find_dfe(R0, p0, c0)` and print result alongside the endemic equilibria

Extend the `@testset "Stability classification"` test to cover the DFE.

---

## Commit 3 — Plot correctness

File: `main.jl`, `src/trajectories.jl`

### B1 — Fix equilibrium marker bug

`main.jl:99`

- Change `Equilibria[1][3]` → `Equilibria[2][3]` in the second `scatter!` call

### B2 — Fix time-series x-axis

`src/trajectories.jl`

- Return `t` from `my_trajectory` alongside `x` and `I` (change return value to `[x, I, t]`)

`main.jl`

- Unpack `t` from `my_trajectory` result
- Pass `t` as first argument to `lines!` in the `ft` time-series plot

---

## Commit 4 — Dead code removal + docs

### B3 — Remove dead duplicate `ft` block

`main.jl:137–159`: delete the second `ft = Figure(...)` block (save already ran at line 133, display is commented out).

### B4 — Remove unused function

`src/trajectories.jl`: delete `my_trajectory_no_condition`.

### B5 — Remove debug `println` calls

- `src/trajectories.jl:11–12`
- `src/equilibria.jl:61`
- `main.jl:77`

### B6 — Remove commented-out exploration code

`main.jl:20–35`: delete the commented-out initial condition grids.

### B7 — Remove obsolete script

```
git rm src/_obsolete_mono_script.jl
```

### A4 — Fix variable description

`README.md`: change "fraction of the population adopting a protective behaviour" to "fraction of the population adopting risky/normal behaviour (non-protective)". Update the parameter table description accordingly.
