# Simple Disease–Behaviour Model

Simulates a coupled SIS epidemic and behavioural dynamics system. Two state variables evolve together:

- **x** — fraction of the population adopting a protective behaviour (e.g. social distancing)
- **I** — fraction of the population infected

Behaviour influences transmission (via `x` in the force of infection), and infection pressure shifts behaviour adoption. The model is analysed for equilibria and simulated as a continuous dynamical system.

## How to run

```julia
julia main.jl
```

Requires Julia with the following packages (installed automatically on first run):
`DynamicalSystems`, `OrdinaryDiffEq`, `GLMakie`, `Roots`

## Output

Plots are saved to `simple-disease-behaviour-plots/`:
- `f1.png` — phase portrait (x vs I)
- `ft.png` — time series of I

## File structure

| File | Purpose |
|---|---|
| `main.jl` | Entry point: loads packages, runs simulation, produces plots |
| `src/params.jl` | Model parameters (R₀, p, c, simulation time, solver settings) |
| `src/model.jl` | ODE right-hand side and initial condition validator |
| `src/equilibria.jl` | Jacobian, eigenvalue analysis, equilibrium finder |
| `src/trajectories.jl` | Trajectory simulation wrappers |

## Key parameters (`params.jl`)

| Parameter | Value | Meaning |
|---|---|---|
| `R0` | 3.0 | Basic reproduction number |
| `p0` | 1.05 | Imitation rate |
| `c0` | 0.19 | Spontaneous behaviour switching rate |
| `T` | 100 | Simulation time |
