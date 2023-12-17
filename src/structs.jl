abstract type AbstractQOEM  <: ContinuousUnivariateDistribution end

"""

    QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM

A model object for a quantum model for medical diagnosis. 
    
# Basis vectors

A person's cognitive state Ψ is represented in a 4D belief space where basis vectors correspond to 
hypotheses and evidence states: 

1. disease present and positive evidence for disease
2. disease present and negative evidence for disease
3. disease absent and positive evidence for disease
4. disease absent and negative evidence for disease

# Fields 

- `Ψ::V`: initial state vector (superposition)
- `γₕ::T`: rotation for positive evidence from medical history
- `γₗ::T`: rotation for negative evidence from laboratory test
- `σ::T`: the standard deviation of probability judgments

# Example 

```julia
Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
γₕ = 4.4045 / √(.5)
γₗ = 0.3306 / √(.5)
σ = .10
model = QOEM(;Ψ, γₕ, γₗ, σ)
predict(model)
```

# References 

Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.
"""
struct QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM
    Ψ::V
    γₕ::T 
    γₗ::T
    σ::T
end

QOEM(;Ψ, γₕ, γₗ, σ) = QOEM(Ψ, γₕ, γₗ, σ)

function QOEM(Ψ, γₕ, γₗ, σ)
    _, γₕ, γₗ, σ = promote(Ψ[1], γₕ, γₗ, σ)
    Ψ = convert(Vector{typeof(γₕ)}, Ψ)
    return QOEM(Ψ, γₕ, γₗ, σ)
end