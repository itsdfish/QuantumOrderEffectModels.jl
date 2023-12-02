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
4. disease present and negative evidence 

# Fields 

- `Ψ::V`: initial state vector (superposition)
- `γₚ::T`: rotation for positive evidence 
- `γₙ::T`: rotation for negative evidence

# Example 

```julia
Ψ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]
γₚ = 4.4045 / √(.5)
γₙ = 0.3306 / √(.5)
dist = QOEM(;Ψ, γₚ, γₙ)
```

# References 

Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.
"""
struct QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM
    Ψ::V
    γₚ::T 
    γₙ::T
end

QOEM(;Ψ, γₚ, γₙ) = QOEM(Ψ, γₚ, γₙ)

function QOEM(Ψ, γₚ, γₙ)
    _, γₚ, γₙ = promote(Ψ[1], γₚ, γₙ)
    Ψ = convert(Vector{typeof(γₚ)}, Ψ)
    return QOEM(Ψ, γₚ, γₙ)
end