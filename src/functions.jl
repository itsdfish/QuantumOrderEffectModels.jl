"""
    predict(model::AbstractQOEM; t = 1)

Returns predicted response probability for the following conditions:

    
# Arguments

- `model::AbstractQOEM`

# Keywords

- `t = 1`: time of decision

# Output 

The output is a vector of predictions corresponding to the following conditions:

1. Pr(disease)
2. Pr(disease | positive)
3. Pr(disease | positive, negative)
4. Pr(disease | negative)
5. Pr(disease | negative, positive)

# Example 

```julia 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₚ = 2
γₙ = .5
σ = .10
model = QOEM(;Ψ, γₚ, γₙ, σ)
predict(model)
```
"""
function predict(model::AbstractQOEM; t = 1.0)
    (;Ψ, γₚ, γₙ) = model

    H1 = make_H1(model)
    H2 = make_H2(model)
    
    H = √(1/2) .* (H1 .+ H2)
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
    Ψn = Uni * Ψ
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
    make_H1(model::AbstractQOEM)

Creates a Hermitian matrix which rotates 

# Arguments 

- `model::AbstractQOEM`: a model object for a quantum order effect model 
"""
make_H1(model::AbstractQOEM) = kron(I(2), [1. 1; 1 -1])


"""
    make_H2(model::AbstractQOEM)

# Arguments 

- `model::AbstractQOEM`: a model object for a quantum order effect model 
"""
function make_H2(model::AbstractQOEM)
    H = kron(fill(1.0, 2, 2), I(2))
    H[2,2] = -1.0
    H[3,3] = -1.0
    return H
end

rand(model::AbstractQOEM; t = 1) = rand(model, 1; t = 1)

"""
    rand(model::AbstractQOEM, n::Int; t = 1)

Generates simulated data for the following conditions:

1. Pr(disease)
2. Pr(disease | positive)
3. Pr(disease | positive, negative)
4. Pr(disease | negative)
5. Pr(disease | negative, positive)

# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 

# Keywords

- `t = 1`: time of decision

# Example 

```julia 
Ψ = @. √([.35,.35,.15,.15])
γₚ = 2
γₙ = .5
σ = .10
model = QOEM(;Ψ, γₚ, γₙ, σ)
data = rand(model, 100)
```
"""
function rand(model::AbstractQOEM, n::Int; t = 1)
    μs = predict(model; t)
    σ = model.σ
    θ = to_beta.(μs, [σ])
    f(Θ) = map(θ -> round_val(rand(Beta(θ...)), .05), Θ)
    return [f(θ) for _ ∈ 1:n]
end

"""
    pdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1)

Returns the joint probability density given data for the following conditions:

1. Pr(disease)
2. Pr(disease | positive)
3. Pr(disease | positive, negative)
4. Pr(disease | negative)
5. Pr(disease | negative, positive)  

# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = 1`: time of decision
"""
function pdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1)
    Θ = predict(model; t)
    return prod(@. pdf(Binomial(n, Θ), n_d)) 
end

"""
    logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1)

Returns the joint log density given data for the following conditions:

1. Pr(disease)
2. Pr(disease | positive)
3. Pr(disease | positive, negative)
4. Pr(disease | negative)
5. Pr(disease | negative, positive)  
    

# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = 1: time of decision

# Example 

```julia 
using QuantumOrderEffectModels 
model = QPDM(;μd=.51, γ=2.09)
n_trials = 100
data = rand(model, n_trials)
logpdf(model, n_trials, data)
```
"""
function logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1)
    Θ = predict(model; t)
    return sum(@. logpdf(Binomial(n, Θ), n_d))
end

loglikelihood(d::AbstractQOEM, data::Tuple) = logpdf(d, data...)

logpdf(model::AbstractQOEM, x::Tuple) = logpdf(model, x...)

function round_val(p, r) 
    v = 1 / r 
    return round(p * v) / v
end

function to_beta(μ, σ)
    α = ((1 - μ) / σ^2 - (1 / μ)) * μ^2
    β = α * ((1 / μ) - 1)
    return α, β
end