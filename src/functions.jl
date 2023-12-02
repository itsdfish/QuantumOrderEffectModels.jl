"""
    predict(dist::AbstractQOEM; t = π / 2)

Returns predicted response probability for the following conditions:
    
# Arguments

- `dist::AbstractQOEM`

# Keywords

- `t = π / 2`: time of decision

# Example 

```julia 
using QuantumPrisonersDilemmaModel 
model = QPDM(;μd=.51, γ=2.09)
predict(model)
```
"""
function predict(dist::AbstractQOEM; t = 1.0)
    (;Ψ, γₚ, γₙ) = dist
    H1 = make_H1(dist)
    H2 = make_H2(dist)
    
    H = H1 .+ H2
    # unitary transformation matrix from initial state to 
    # positive evidence
    Upi = exp(-im * t * γₚ * H)
    # unitary transformation matrix from initial state to 
    # positive evidence
    Uni = exp(-im * t * γₙ * H)
    
    #  projector for positive evidence
    Pp = Diagonal([1,0,1,0])
    # projector for positive evidence and disease present
    Pd = Diagonal([1,1,0,0])
    #  projector for negative evidence
    Pn = Diagonal([0,1,0,1])
    
    # probability of disease before observing additional evidence
    proj_d = Pd * Ψ
    prob_d = proj_d' * proj_d
    
    # initial state in terms of positive evidence
    Ψp = Upi * Ψ
    # state projected onto observed positive evidence
    Ψp′ = Pp * Ψp
    # update/normalize the state
    Ψp′ = Ψp′ ./ norm(Ψp′)
    # projection onto disease present basis
    proj_d = Pd * Ψp′
    # probability of disease given positive evidence 
    prob_dgp = real(proj_d'* proj_d)
    
    # unitary transformsations: initial to positive, then to negative
    Ψpn = Uni * Upi'* Ψp′
    Ψpn′ = Pn * Ψpn
    Ψpn′ = Ψpn′ ./ norm(Ψpn′)
    # projection onto disease present basis
    proj_d = Pd * Ψpn′
    # probability of disease given positive evidence then negative evidence
    prob_dgpn = real(proj_d'* proj_d)
    
    # initial state in terms of negative evidence
    Ψp = Uni * Ψ
    # state projected onto observed negative evidence
    Ψn′ = Pn * Ψn
    # update/normalize the state
    Ψn′ = Ψn′ ./ norm(Ψn′)
    # projection onto disease present basis
    proj_d = Pd * Ψn′
    # probability of disease given negative evidence 
    prob_dgn = real(proj_d'* proj_d)
    
    # unitary transformsations: initial to negative, then to positive
    Ψpp = Upi * Uni'* Ψn′
    Ψpp′ = Pp * Ψpp
    Ψpp′ = Ψpp′ ./ norm(Ψpp′)
    # projection onto disease present basis
    proj_d = Pd * Ψpp′
    # probability of disease given negative and then positive evidence 
    prob_dgnp = real(proj_d'* proj_d)
    
    return [prob_d, prob_dgp, prob_dgpn, prob_dgn, prob_dgnp]    
end

"""
    make_H1(dist::AbstractQOEM)

Creates a Hermitian matrix which rotates in favor of defecting or cooperating depending on 
μd and μd. 

# Arguments 

- `μd`: utility for defecting 
- `μc`: utility for cooperating
"""
make_H1(dist::AbstractQOEM) = kron(I(2), [1. 1; 1 -1])


"""
    make_H2(γ)

# Arguments 

- `γ`: entanglement parameter which aligns beliefs and actions
"""
function make_H2(dist::AbstractQOEM)
    H = kron(fill(1.0, 2, 2), I(2))
    H[2,2] = -1.0
    H[3,3] = -1.0
    return H
end

rand(dist::AbstractQOEM; t = π / 2) = rand(dist, 1; t = π / 2)

"""
    rand(dist::AbstractQOEM, n::Int; t = π / 2)

Generates simulated data for the following conditions:

1. Player 2 is told that player 1 defected
2. Player 2 is told that player 1 cooperated
3. Player 2 is not informed of player 1's action

# Arguments

- `dist::AbstractQOEM`
- `n`: the number of trials per condition 

# Keywords

- `t = π / 2`: time of decision

# Example 

```julia 
using QuantumPrisonersDilemmaModel 
model = QPDM(;μd=.51, γ=2.09)
data = rand(model, 100)
```
"""
function rand(dist::AbstractQOEM, n::Int; t = π / 2)
    Θ = predict(dist; t)
    return @. rand(Binomial(n, Θ))
end

"""
    pdf(dist::AbstractQOEM, n::Int, n_d::Vector{Int}; t = π / 2)

Returns the joint probability density given data for the following conditions:

1. Player 2 is told that player 1 defected
2. Player 2 is told that player 1 cooperated
3. Player 2 is not informed of player 1's action
    

# Arguments

- `dist::AbstractQOEM`
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = π / 2`: time of decision
"""
function pdf(dist::AbstractQOEM, n::Int, n_d::Vector{Int}; t = π / 2)
    Θ = predict(dist; t)
    return prod(@. pdf(Binomial(n, Θ), n_d)) 
end

"""
    logpdf(dist::AbstractQOEM, n::Int, n_d::Vector{Int}; t = π / 2)

Returns the joint log density given data for the following conditions:

1. Player 2 is told that player 1 defected
2. Player 2 is told that player 1 cooperated
3. Player 2 is not informed of player 1's action

# Arguments

- `dist::AbstractQOEM`
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = π / 2`: time of decision

# Example 

```julia 
using QuantumPrisonersDilemmaModel 
model = QPDM(;μd=.51, γ=2.09)
n_trials = 100
data = rand(model, n_trials)
logpdf(model, n_trials, data)
```
"""
function logpdf(dist::AbstractQOEM, n::Int, n_d::Vector{Int}; t = π / 2)
    Θ = predict(dist; t)
    return sum(@. logpdf(Binomial(n, Θ), n_d))
end

loglikelihood(d::AbstractQOEM, data::Tuple) = logpdf(d, data...)

logpdf(dist::AbstractQOEM, x::Tuple) = logpdf(dist, x...)