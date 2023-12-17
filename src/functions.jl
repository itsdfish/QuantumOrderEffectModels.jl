"""
    predict(model::AbstractQOEM; t = 1)

Returns predicted response probability for the following conditions:
    
# Arguments

- `model::AbstractQOEM`

# Keywords

- `t = 1`: time of decision

# Output 

The output is a vector of predictions corresponding to the following conditions:

## Order 1
1. Pr(disease)
2. Pr(disease | lab test)
3. Pr(disease | lab test, medical history)

## Order 2
4. Pr(disease)
5. Pr(disease | medical history)
6. Pr(disease | medical history, lab test)

# Example 

```julia 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₕ = 2
γₗ = .5
σ = .05
model = QOEM(;Ψ, γₕ, γₗ, σ)
predict(model)
```
"""
function predict(model::AbstractQOEM; t = 1.0)
    (;Ψ, γₕ, γₗ) = model

    H1 = make_H1(model)
    H2 = make_H2(model)
    H = √(1/2) .* (H1 .+ H2)

    # unitary transformation matrix for positive evidence from medical history
    Upi = exp(-im * t * γₕ * H)
    # unitary transformation matrix for negative evidence from laboratory test
    Uni = exp(-im * t * γₗ * H)
    
    #  projector for positive evidence
    Pp = Diagonal([1,0,1,0])
    # projector for disease present
    Pd = Diagonal([1,1,0,0])
    # projector for negative evidence
    Pn = Diagonal([0,1,0,1])
    
    # probability of disease before observing additional evidence
    proj_d = Pd * Ψ
    prob_d = proj_d' * proj_d
    
    # state projected onto observed positive evidence from medical history
    Ψp = Pp * Upi * Ψ
    # update/normalize the state
    Ψp = Ψp ./ norm(Ψp)
    # projection onto disease present basis
    proj_d = Pd * Ψp
    # probability of disease given positive evidence 
    prob_dgp = real(proj_d'* proj_d)
    
    # project onto negative evidence basis for laboratory test
    Ψpn = Pn * Uni * Upi'* Ψp
    # update/normalize the state
    Ψpn = Ψpn ./ norm(Ψpn)
    # projection onto disease present basis
    proj_d = Pd * Ψpn
    # probability of disease given positive evidence then negative evidence
    prob_dgpn = real(proj_d' * proj_d)
    
    # state projected onto observed negative evidence for laboratory test
    Ψn = Pn * Uni * Ψ
    # update/normalize the state
    Ψn = Ψn ./ norm(Ψn)
    # projection onto disease present basis
    proj_d = Pd * Ψn
    # probability of disease given negative evidence 
    prob_dgn = real(proj_d' * proj_d)
    
    # unitary transformsations: negative to initial, then to positive
    Ψnp = Pp * Upi * Uni'* Ψn
    # update/normalize the state
    Ψnp = Ψnp ./ norm(Ψnp)
    # projection onto disease present basis
    proj_d = Pd * Ψnp
    # probability of disease given negative and then positive evidence 
    prob_dgnp = real(proj_d' * proj_d)
    
    return [prob_d, prob_dgp, prob_dgpn, prob_d, prob_dgn, prob_dgnp]    
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

rand(model::AbstractQOEM; t = 1, r = .05) = rand(model, 1; t, r)

"""
    rand(model::AbstractQOEM, n::Int; t = 1, r = .05)

Generates simulated data for the following conditions:

## Order 1
1. Pr(disease)
2. Pr(disease | lab test)
3. Pr(disease | lab test, medical history)

## Order 2
4. Pr(disease)
5. Pr(disease | medical history)
6. Pr(disease | medical history, lab test)

# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 

# Keywords

- `t = 1`: time of decision

# Example 

```julia 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₕ = 2
γₗ = .5
σ = .05
n_trials = 100
model = QOEM(;Ψ, γₕ, γₗ, σ)
data = rand(model, n_trials)
```
"""
function rand(model::AbstractQOEM, n::Int; t = 1, r = .05)
    μs = predict(model; t)
    σ = model.σ
    θ = to_beta.(μs, [σ])
    #f(Θ) = map(θ -> round_val(rand(Beta(θ...)), r), Θ)
    f(Θ) = map(θ -> rand(Beta(θ...)), Θ)
    return [f(θ) for _ ∈ 1:n]
end

"""
    logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1, r = .05)

Returns the joint log density given data for the following conditions:

## Order 1
1. Pr(disease)
2. Pr(disease | lab test)
3. Pr(disease | lab test, medical history)

## Order 2
4. Pr(disease)
5. Pr(disease | medical history)
6. Pr(disease | medical history, lab test)
    
# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = 1: time of decision

# Example 

```julia 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₕ = 2
γₗ = .5
σ = .05
n_trials = 100
model = QOEM(;Ψ, γₕ, γₗ, σ)
data = rand(model, n_trials)
logpdf(model, data)
```
"""
function pdf(model::AbstractQOEM, data::Vector{Vector{T}}; t = 1, r = .05) where {T<:Real}
    Θ = predict(model; t)
    L = 1.0
    for d ∈ data 
        L *= _pdf(model, Θ, d; r)
    end
    return L
end

function _pdf(model::AbstractQOEM, Θ, data; r = .05)
    L = 1.0
    σ = model.σ
    for i ∈ 1:length(Θ) 
        parms = to_beta(Θ[i], σ)
        lb,ub = get_bounds(Θ[i], r)
        L *= resp_prob(parms, lb, ub)
    end
    return L
end

"""
    logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1, r = .05)

Returns the joint log density given data for the following conditions:

## Order 1
1. Pr(disease)
2. Pr(disease | lab test)
3. Pr(disease | lab test, medical history)

## Order 2
4. Pr(disease)
5. Pr(disease | medical history)
6. Pr(disease | medical history, lab test)
    
# Arguments

- `model::AbstractQOEM`: a model object for a quantum order effect model 
- `n`: the number of trials per condition 
- `n_d`: the number of defections in each condition 

# Keywords

- `t = 1: time of decision

# Example 

```julia 
Ψ = @. √([.35,.35,.15,.15])
γₕ = 2
γₗ = .5
σ = .05
n_trials = 100
model = QOEM(;Ψ, γₕ, γₗ, σ)
data = rand(model, n_trials)
logpdf(model, data)
```
"""
function logpdf(model::AbstractQOEM, data::Vector{Vector{T}}; t = 1, r = .05) where {T<:Real}
    Θ = predict(model; t)
    LL = 0.0
    for d ∈ data 
        LL += _logpdf(model, Θ, d; r)
    end
    return LL
end

function _logpdf(model::AbstractQOEM, Θ, data; r = .05)
    LL = 0.0
    σ = model.σ
    for i ∈ 1:length(Θ) 
        parms = to_beta(Θ[i], σ)
        lb,ub = get_bounds(Θ[i], r)
        #LL += log(resp_prob(parms, lb, ub))
        LL += logpdf(Beta(parms...), data[i])
    end
    return LL
end

loglikelihood(d::AbstractQOEM, data::Vector{Vector{T}}) where {T<:Real} = logpdf(d, data)

function round_val(p, r) 
    v = 1 / r 
    return round(p * v) / v
end

function to_beta(μ, σ)
    α = ((1 - μ) / σ^2 - (1 / μ)) * μ^2
    β = α * ((1 / μ) - 1)
    α = max(α, eps())
    β = max(β, eps())
    return α, β
end

get_bounds(x, r) = max(x - r / 2, 0.0), min(x + r / 2, 1.0)

function resp_prob(parms, lb, ub) 
    dist = Beta(parms...)
    return cdf(dist, ub) - cdf(dist, lb)
end

function mylogpdf(α, β, x)
    lb,ub = get_bounds(x, .05)
    return log(resp_prob([α,β], lb, ub))
end