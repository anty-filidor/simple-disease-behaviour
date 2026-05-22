using Test
using LinearAlgebra
using OrdinaryDiffEq
using Roots

include("../src/params.jl")
include("../src/model.jl")
include("../src/equilibria.jl")

@testset "ODE RHS" begin
    R, p, c = R0, p0, c0
    xs = calculate_values(R, p, c)
    @assert length(xs) >= 1 "No equilibria found for default params"
    x_star = xs[1]
    I_star = (R * x_star - 1) / (R * x_star)

    dxdt =
        -p*(1+I_star)*x_star*(1-x_star)^2 +
        p*(1-I_star)*(1-x_star)*x_star^2 +
        c*(1-I_star)*(1-x_star)
    dIdt = R * x_star * I_star * (1 - I_star) - I_star

    @test dxdt ≈ 0.0 atol=1e-10
    @test dIdt ≈ 0.0 atol=1e-10

    I_perturbed = I_star + 0.1
    dIdt_perturbed = R * x_star * I_perturbed * (1 - I_perturbed) - I_perturbed
    @test dIdt_perturbed < 0
end

@testset "Boundary validator" begin
    @test cond(0.5, 0.3) == true
    @test cond(-0.1, 0.3) == false
    @test cond(0.5, 0.0) == false
    @test cond(0.5, 1.1) == false
end

@testset "Equilibrium finder" begin
    R, p, c = R0, p0, c0
    roots = calculate_values(R, p, c)
    @test length(roots) >= 1
    @test all(0 .< roots .< 1)

    for x_star in roots
        I_star = (R * x_star - 1) / (R * x_star)
        dxdt =
            -p*(1+I_star)*x_star*(1-x_star)^2 +
            p*(1-I_star)*(1-x_star)*x_star^2 +
            c*(1-I_star)*(1-x_star)
        dIdt = R * x_star * I_star * (1 - I_star) - I_star
        @test dxdt ≈ 0.0 atol=1e-10
        @test dIdt ≈ 0.0 atol=1e-10
    end

    # No valid endemic equilibria when R is below the invasion threshold
    @test length(find_equilibria(1.0, p, c)) == 0
end

@testset "Jacobian correctness" begin
    x_t, z_t = 0.4, 0.2
    R, p, c = R0, p0, c0

    J = jacob(x_t, z_t, R, p, c)

    @test J[1, 1] ≈ -p*(1 - 6x_t + 6x_t^2 + z_t - 2x_t*z_t) - c*(1-z_t) atol=1e-10
    @test J[1, 2] ≈ -(1-x_t)*(p*x_t + c) atol=1e-10
    @test J[2, 1] ≈ R*z_t*(1-z_t) atol=1e-10
    @test J[2, 2] ≈ R*x_t*(1-2z_t) - 1 atol=1e-10
end

@testset "Stability classification" begin
    R, p, c = R0, p0, c0
    equilibria = find_equilibria(R, p, c)
    @test length(equilibria) > 0

    for eq in equilibria
        x_star, I_star, is_stable = eq
        vals, _ = compute_eigenvalues_jacob(x_star, I_star, R, p, c)
        if is_stable
            @test all(real.(vals) .<= 0)
        end
    end

    _, _, dfe_stable = find_dfe(R, p, c)
    @test dfe_stable == false

    _, _, dfe_stable_low = find_dfe(0.5, p, c)
    @test dfe_stable_low == true
end
