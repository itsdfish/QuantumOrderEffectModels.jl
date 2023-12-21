using SafeTestsets

@safetestset "predict" begin
    using QuantumOrderEffectModels
    using Test

    Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
    γₕ = 4.4045 / √(.5)
    γₗ = 0.3306 / √(.5)
    σ = .10
    dist = QOEM(;Ψ, γₕ, γₗ, σ)
    preds = predict(dist)
    true_preds = [0.676, 0.793, 0.504, 0.676, 0.437, 0.59]

    @test preds ≈ true_preds atol = 1e-3
end

@safetestset "to_beta" begin
    using Distributions
    using QuantumOrderEffectModels
    using QuantumOrderEffectModels: to_beta
    using Test
    
    α = 2
    β = 3 
    dist = Beta(α, β)
    μ = mean(dist)
    σ = std(dist)
    α′,β′ = to_beta(μ, σ)

    @test α ≈ α′
    @test β ≈ β′
end

@safetestset "rand" begin
    @safetestset "rand 1" begin
        using QuantumOrderEffectModels
        using Test
        using Random 
        using Statistics

        Random.seed!(7878)
    
        Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
        γₕ = 4.4045 / √(.5)
        γₗ = 0.3306 / √(.5)
        σ = .10
        dist = QOEM(;Ψ, γₕ, γₗ, σ)
        preds = predict(dist)
        data = rand(dist, 10_000)

        @test preds ≈ mean(data) rtol = .01
        @test fill(σ, 6) ≈ std(data) rtol = .01
    end

    @safetestset "rand 2" begin
        using QuantumOrderEffectModels
        using Test
        using Random 
        using Statistics

        Random.seed!(52)
    
        Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
        γₕ = 3.5
        γₗ = -2.0
        σ = .15
        dist = QOEM(;Ψ, γₕ, γₗ, σ)
        preds = predict(dist)
        data = rand(dist, 10_000)

        @test preds ≈ mean(data) rtol = .01
        @test fill(σ, 6) ≈ std(data) rtol = .01
    end
end

@safetestset "get_bounds" begin
    using QuantumOrderEffectModels: get_bounds
    using Test

    lb,ub = get_bounds(.4, .05)
    @test lb ≈ .375
    @test ub ≈ .425

    lb,ub = get_bounds(.0, .05)
    @test lb ≈ .0
    @test ub ≈ .025

    lb,ub = get_bounds(1, .05)
    @test lb ≈ .975
    @test ub ≈ 1
end

@safetestset "logpdf" begin
    using QuantumOrderEffectModels
    using Random
    using Test

    Random.seed!(584)

    Θ = (Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)],
        γₕ = 4.0,
        γₗ = 2.0, 
        σ = .05)
    dist = QOEM(;Θ...)
    r = .05
    data = rand(dist, 10_000; r)

    γₕs = range(Θ.γₕ * .80, Θ.γₕ * 1.20, length = 100)
    LLs = map(γₕ -> logpdf(QOEM(;Θ..., γₕ), data; r), γₕs)
    _,max_idx = findmax(LLs)
    @test γₕs[max_idx] ≈ Θ.γₕ rtol = .01

    γₗs = range(Θ.γₗ * .80, Θ.γₗ * 1.20, length = 100)
    LLs = map(γₗ -> logpdf(QOEM(;Θ..., γₗ), data), γₗs)
    _,max_idx = findmax(LLs)
    @test γₗs[max_idx] ≈ Θ.γₗ rtol = .01
end