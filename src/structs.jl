abstract type AbstractQOEM  <: ContinuousUnivariateDistribution end

"""

    QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM

A model object for a quantum model for medical diagnosis. 
    
# Basis vectors

A person's cognitive state Ψ is represented in a 4D belief space where basis vectors correspond to 
hypotheses and evidence states: 

1. disease present and positive evidence 
2. disease present and negative evidence 
3. disease absent and positive evidence 
4. disease absent and negative evidence 

# Fields 

- `Ψ::V`: initial state vector (superposition)
- `γₚ::T`: rotation for positive evidence 
- `γₙ::T`: rotation for negative evidence
- `σ::T`: the standard deviation of probability judgments

# Example 

```julia
Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
γₚ = 4.4045 / √(.5)
γₙ = 0.3306 / √(.5)
σ = .10
model = QOEM(;Ψ, γₚ, γₙ, σ)
predict(model)
```

# References 

Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.
"""
struct QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM
    Ψ::V
    γₚ::T 
    γₙ::T
    σ::T
end

QOEM(;Ψ, γₚ, γₙ, σ) = QOEM(Ψ, γₚ, γₙ, σ)

function QOEM(Ψ, γₚ, γₙ, σ)
    _, γₚ, γₙ, σ = promote(Ψ[1], γₚ, γₙ, σ)
    Ψ = convert(Vector{typeof(γₚ)}, Ψ)
    return QOEM(Ψ, γₚ, γₙ, σ)
end