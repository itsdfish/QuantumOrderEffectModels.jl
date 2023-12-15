# Overview

This page provides an overview of the API along with examples. 

# Make Predictions

The quantum prisoner's dilemma model (QPDM) generates predictions for three conditions:

1. Player 1 is told that player 2 defected
2. Player 1 is told that player 2 cooperated
3. Player 1 is not informed of the action of player 2

```@example 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₚ = 2
γₙ = .5
σ = .05
model = QOEM(;Ψ, γₚ, γₙ, σ)
predict(model)
```

# Simulate Model

The code block below demonstrates how to generate simulated data from the model using `rand`. In the example, we will generate 100 simulated trials for each condition. 
```@example 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₚ = 2
γₙ = .5
σ = .05
n_trials = 100
model = QOEM(;Ψ, γₚ, γₙ, σ)
data = rand(model, n_trials)
```

# Evaluate Log Likelihood

The log likelihood of data can be evaluated using `logpdf`. In the code block below, we generate simulated data and evaluate the logpdf: 
```@example 
using QuantumOrderEffectModels
Ψ = @. √([.35,.35,.15,.15])
γₚ = 2
γₙ = .5
σ = .05
n_trials = 100
model = QOEM(;Ψ, γₚ, γₙ, σ)
data = rand(model, n_trials)
logpdf(model, data)
```