abstract type AbstractQOEM  <: ContinuousUnivariateDistribution end

"""

    QOEM{T<:Real} <: AbstractQOEM

A model object for the Quantum Prisoner's Dilemma Model. The QOEM has four basis states:
    
1. opponent defects and you defect 
2. opponent defects and you cooperate 
3. opponent cooperates and you defect 
4. opponent cooperates and you cooperate

The bases are orthonormal and in standard form. The model assumes three conditions:

1. Player 2 is told that player 1 defected
2. Player 2 is told that player 1 cooperated
3. Player 2 is not informed of the action of player 1


Model inputs and outputs are assumed to be in the order above. 

# Fields 

- `μd`: utility for defecting 
- `μc`: utility for cooperating 
- `γ`: entanglement parameter for beliefs and actions 

# Example 

```julia
using QuantumPrisonersDilemmaModel
model = QOEM(;μd=.51, γ=2.09)
```

# References 

Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.
"""
struct QOEM{T<:Real} <: AbstractQOEM
    Ψ::AbstractVector{<:T}
    γₚ::T 
    γₙ::T
end

QOEM(;Ψ, γₚ, γₙ) = QOEM(Ψ, γₚ, γₙ)

function QOEM(Ψ, γₚ, γₙ)
    _, γₚ, γₙ = promote(Ψ[1], γₚ, γₙ)
    Ψ = convert(Vector{typeof(γₚ)}, Ψ)
    return QOEM(Ψ, γₚ, γₙ)
end