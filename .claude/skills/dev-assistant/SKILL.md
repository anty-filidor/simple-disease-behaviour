---
name: dev-assistant
description: >
  Development assistant for this Julia disease-behaviour modelling repository. Use this skill whenever
  working on any development task here — writing code, debugging, running simulations, adding features,
  managing dependencies, or committing changes. Provides coding conventions, active workflow guidance,
  and project standards.
---

## Session start

At the start of every session, read `README.md` in the working directory root to get repo-specific
context: purpose, setup instructions, usage, and available commands.

---

## Language

Use **British English** in all text: comments, docstrings, commit messages, and documentation.

---

## Code style

- Line length: **100 characters**
- **Type annotations**: use Julia type annotations on function signatures where they clarify intent;
  follow the pattern already established in the file
- **Naming**: `snake_case` for functions and variables, `PascalCase` for types and modules,
  `UPPER_SNAKE_CASE` for module-level constants
- **Docstrings**: `"""..."""` above the function definition, Julia stdlib style:
  - First line: indented call signature with return type (`    f(x, R0) -> Float64`)
  - Blank line, then a single concise sentence describing what the function does
  - If parameters need clarification, add a brief `# Arguments` section — otherwise omit it
  - One-liner variant (no signature line needed): just the sentence on one line

**Docstring examples:**

```julia
"""Short description on one line."""
function retry_count(attempts::Int)

"""
    find_equilibria(R0, p0, c0) -> Vector

Return all equilibria of the system for the given parameters.
"""
function find_equilibria(R0, p0, c0)

"""
    my_trajectory(x, I, R0, p0, c0) -> Union{Tuple, Nothing}

Simulate a trajectory from initial state `(x, I)` and return `(xs, Is)`, or `nothing` if the
initial condition is invalid.

# Arguments
- `x`: initial behavioural adoption fraction ∈ [0, 1]
- `I`: initial infected fraction ∈ [0, 1]
"""
function my_trajectory(x, I, R0, p0, c0)
```

---

## Available CLI tools

| Tool     | Purpose                                                              |
| -------- | -------------------------------------------------------------------- |
| `julia`  | Run scripts and manage packages (`julia --project=. main.jl`)        |
| `git`    | Version control                                                      |

---

## Package management

Add a dependency from the shell:

```bash
julia --project=. -e 'using Pkg; Pkg.add("PackageName")'
```

Never edit `Project.toml` or `Manifest.toml` by hand for dependency changes.

---

## Running the project

```bash
julia --project=. main.jl
```

Instantiate the environment on a fresh clone:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## New branch workflow

Before making any other changes on a new branch:

1. Check the current `version` in `Project.toml` on the parent branch
2. Compare it to `version` on the current branch — if they differ, the bump has already been done;
   skip to step 3. If they are the same, bump `version` in `Project.toml` (patch → minor → major,
   following semver)
3. Then proceed with the branch work

---

## Git workflow

- Move files with `git mv`, never bare `mv`
- Do not commit known-broken code — confirm `julia --project=. main.jl` completes without error and
  pre-commit checks pass
- Commit messages: short imperative subject line in British English
- No co-authorship trailers in commit messages

---

## Pre-commit

Run before every commit:

```bash
pre-commit run --files <changed files>
```

Hooks configured: trailing whitespace, end-of-file fixer, YAML check, case-conflict check, large-file check, and JuliaFormatter (auto-formats `.jl` files in place).

---

## Security

- Do not hardcode secrets, tokens, or passwords in source files — use environment variables
- `rm -rf` and `git push --force` / `git push -f` are blocked at the harness level
