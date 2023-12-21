var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumOrderEffectModels]\nOrder   = [:type, :function]\nPrivate = false","category":"page"},{"location":"api/#QuantumOrderEffectModels.QOEM","page":"API","title":"QuantumOrderEffectModels.QOEM","text":"QOEM{T<:Real,V<:AbstractVector{T}} <: AbstractQOEM\n\nA model object for a quantum model for medical diagnosis. \n\nBasis vectors\n\nA person's cognitive state Ψ is represented in a 4D belief space where basis vectors correspond to  hypotheses and evidence states: \n\ndisease present and positive evidence for disease\ndisease present and negative evidence for disease\ndisease absent and positive evidence for disease\ndisease absent and negative evidence for disease\n\nFields\n\nΨ::V: initial state vector (superposition)\nγₕ::T: rotation for positive evidence from medical history\nγₗ::T: rotation for negative evidence from laboratory test\nσ::T: the standard deviation of probability judgments\n\nExample\n\nΨ = [√(.676 / 2),√(.676 / 2),√(.324 / 2),√(.324 / 2)]\nγₕ = 4.4045 / √(.5)\nγₗ = 0.3306 / √(.5)\nσ = .10\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\npredict(model)\n\nReferences\n\nTrueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.rand-Tuple{QuantumOrderEffectModels.AbstractQOEM, Int64}","page":"API","title":"Base.rand","text":"rand(model::AbstractQOEM, n::Int; t = 1, r = .05)\n\nGenerates simulated data for the following conditions:\n\nOrder 1\n\nPr(disease)\nPr(disease | lab test)\nPr(disease | lab test, medical history)\n\nOrder 2\n\nPr(disease)\nPr(disease | medical history)\nPr(disease | medical history, lab test)\n\nArguments\n\nmodel::AbstractQOEM: a model object for a quantum order effect model \nn: the number of trials per condition \n\nKeywords\n\nt = 1: time of decision\n\nExample\n\nusing QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nn_trials = 100\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\ndata = rand(model, n_trials)\n\n\n\n\n\n","category":"method"},{"location":"api/#Distributions.logpdf-Union{Tuple{T}, Tuple{QuantumOrderEffectModels.AbstractQOEM, Array{Vector{T}, 1}}} where T<:Real","page":"API","title":"Distributions.logpdf","text":"logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1, r = .05)\n\nReturns the joint log density given data for the following conditions:\n\nOrder 1\n\nPr(disease)\nPr(disease | lab test)\nPr(disease | lab test, medical history)\n\nOrder 2\n\nPr(disease)\nPr(disease | medical history)\nPr(disease | medical history, lab test)\n\nArguments\n\nmodel::AbstractQOEM: a model object for a quantum order effect model \nn: the number of trials per condition \nn_d: the number of defections in each condition \n\nKeywords\n\n`t = 1: time of decision\n\nExample\n\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nn_trials = 100\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\ndata = rand(model, n_trials)\nlogpdf(model, data)\n\n\n\n\n\n","category":"method"},{"location":"api/#Distributions.pdf-Union{Tuple{T}, Tuple{QuantumOrderEffectModels.AbstractQOEM, Array{Vector{T}, 1}}} where T<:Real","page":"API","title":"Distributions.pdf","text":"logpdf(model::AbstractQOEM, n::Int, n_d::Vector{Int}; t = 1, r = .05)\n\nReturns the joint log density given data for the following conditions:\n\nOrder 1\n\nPr(disease)\nPr(disease | lab test)\nPr(disease | lab test, medical history)\n\nOrder 2\n\nPr(disease)\nPr(disease | medical history)\nPr(disease | medical history, lab test)\n\nArguments\n\nmodel::AbstractQOEM: a model object for a quantum order effect model \nn: the number of trials per condition \nn_d: the number of defections in each condition \n\nKeywords\n\n`t = 1: time of decision\n\nExample\n\nusing QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nn_trials = 100\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\ndata = rand(model, n_trials)\nlogpdf(model, data)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumOrderEffectModels.predict-Tuple{QuantumOrderEffectModels.AbstractQOEM}","page":"API","title":"QuantumOrderEffectModels.predict","text":"predict(model::AbstractQOEM; t = 1)\n\nReturns predicted response probability for the following conditions:\n\nArguments\n\nmodel::AbstractQOEM\n\nKeywords\n\nt = 1: time of decision\n\nOutput\n\nThe output is a vector of predictions corresponding to the following conditions:\n\nOrder 1\n\nPr(disease)\nPr(disease | lab test)\nPr(disease | lab test, medical history)\n\nOrder 2\n\nPr(disease)\nPr(disease | medical history)\nPr(disease | medical history, lab test)\n\nExample\n\nusing QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\npredict(model)\n\n\n\n\n\n","category":"method"},{"location":"parameter_estimation/#Parameter-Estimation","page":"Parameter Estimation","title":"Parameter Estimation","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"This brief tutorial explains how to performance Bayesian parameter estimation of the QPDM using Pigeons.jl. One complication in estimating the parameters of the QPDM is that the posterior distributions may have multiple modes, which leads to convergence problems with most MCMC algorithms. Pigeons.jl uses a special type of parallel tempering to overcome this challenge. An additional advantage of using Pigeons.jl is the ability to compute Bayes factors from the log marginal likelihood using the function stepping_stone.","category":"page"},{"location":"parameter_estimation/#Load-Packages","page":"Parameter Estimation","title":"Load Packages","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"First, we will load the required packages below. ","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"using Pigeons\nusing QuantumOrderEffectModels\nusing Random\nusing StatsPlots\nusing Turing","category":"page"},{"location":"parameter_estimation/#Generate-Simulated-Data","page":"Parameter Estimation","title":"Generate Simulated Data","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"The next step is to generate some simulated data from which the parameters can be estimated. In the code block below, the utility parameter mu_d is set to one and the entanglement parameter is set to gamma = 2.  A total of 50 trials is generated for each of the three conditions. The resulting values represent the number of defections per condition out of 50.","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Random.seed!(16)\nΨ = @. √([.35,.35,.15,.15])\nparms = (Ψ, γₕ = 2.0, γₗ = 1.0, σ = .05)\nn_trials = 50\nmodel = QOEM(;parms...)\ndata = rand(model, n_trials)","category":"page"},{"location":"parameter_estimation/#Define-Turing-Model","page":"Parameter Estimation","title":"Define Turing Model","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"The next step is to define a Turing model with the @model macro. We will estimate the entanglement parameters using the prior gamma_j sim mathrmnormal(03). The other parameters will be fixed to the data generating values defined in the code block above.","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"@model function turing_model(data, parms)\n    γₕ ~ Normal(0, 3)\n    γₗ ~ Normal(0, 3)\n    σ ~ LogNormal(-1, 1)\n    data ~ QOEM(;parms..., γₕ, γₗ, σ)\nend","category":"page"},{"location":"parameter_estimation/#Estimate-Parameters","page":"Parameter Estimation","title":"Estimate Parameters","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"To estimate the parameters, we need to pass the Turing model to pigeons. The second command converts the output to an MCMCChain object, which can be used for plotting","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"pt = pigeons(\n    target=TuringLogPotential(turing_model(data, parms)), \n    record=[traces],\n    multithreaded=true)\nsamples = Chains(sample_array(pt), [\"γₕ\", \"γₗ\",\"σ\"])\nplot(samples)","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"The trace of the pigeon's sampler is given below:","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"────────────────────────────────────────────────────────────────────────────\n  scans        Λ      log(Z₁/Z₀)   min(α)     mean(α)    min(αₑ)   mean(αₑ) \n────────── ────────── ────────── ────────── ────────── ────────── ──────────\n        2        4.9       -150          0      0.456      0.857      0.961 \n        4       2.47        402   1.61e-23      0.726      0.973      0.997 \n        8       5.03        444   0.000401      0.441          1          1 \n       16       5.27        459      0.136      0.415      0.993      0.999 \n       32       6.11        461     0.0646      0.321          1          1 \n       64       5.51        466      0.232      0.388          1          1 \n      128       5.13        470      0.257       0.43      0.997          1 \n      256       5.18        469      0.311      0.425      0.998          1 \n      512       5.14        469      0.347      0.429      0.999          1 \n 1.02e+03       5.14        469      0.398      0.429      0.998          1 \n────────────────────────────────────────────────────────────────────────────","category":"page"},{"location":"parameter_estimation/#Plot-Posterior-Distribution","page":"Parameter Estimation","title":"Plot Posterior Distribution","text":"","category":"section"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Now we can plot the posterior distribution of gamma with plot. The posterior distribution of gamma has a primary mode around 1 and secondary modes around 2 and 3.5.","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"plot(samples)","category":"page"},{"location":"parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"(Image: )","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"warning: Warning\nThis page is under construction.","category":"page"},{"location":"model_description/#Introduction","page":"Model Description","title":"Introduction","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Prior research has found that the judged likelihood of a hypothesis often depends on the order in which evidence is presented. In other words, the final judgment that a hypothesis is true is different for evidence sequence S_i S_j and evidence sequence S_j S_i. The goal of this tutorial is to describe a quantum order effect model (QOEM) as it applies to a medical diagnosis task. The basic model can be adapted to other tasks in which evidence is evaluated sequentially.  ","category":"page"},{"location":"model_description/#Medical-Diagnosis-Task","page":"Model Description","title":"Medical Diagnosis Task","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Medical doctors routinely perform medical diagnosis as part of their job requirements. Medical diagnosis is the process of assessing the probability a patient has a disease based on symptoms and medical tests. Below, we use the QOEM to understand a how a person evaluates the hypothesis that a fictitous patient has a disease (d) after evidence from two informatino sources are presented sequentially. The judgments are made in two conditions, which vary the presentation order of the evidence. In one condition, a person makes the following judgments:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"the probability of disease based on initial symptoms \nupdate the probability of disease after the medical history S_i provides positive evidence for the disease\nupdate the probability of disease after the laboratory test S_j provides negative evidence for the disease","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The procedure for the other condition is identical except the order of the information sources S_iS_j is reversed: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"the probability of disease based on initial symptoms \nupdate the probability of disease after the laboratory test S_j provides negative evidence for the disease\nupdate the probability of disease after the medical history S_i provides positive evidence for the disease","category":"page"},{"location":"model_description/#Order-Effect","page":"Model Description","title":"Order Effect","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"An order effect occurs when the final probability judgment of disease depends on the order in which evidence is presented. An order effect for two sources of evidence can be stated formally as: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Pr(D=d mid S_i=s_i S_j=s_j) ne Pr(D=d mid S_j=s_j S_i=s_i)","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"where s_is_j in 1-1 and 1 denotes positive evidence for the disease and -1 denotes negative evidence for the disease. Order effects are challenging to explain with classical probability models because evidence is assumed to be commutative (i.e., invariant to order). One workaround for classical probability is to augment the sample space with an additional event representing evidence order. However, augmenting the sample space comes at the cost of creating additional parameters which need to be specified and justified, or estimated from limited data.","category":"page"},{"location":"model_description/#Quantum-Order-Effect-Model","page":"Model Description","title":"Quantum Order Effect Model","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"According to the QOEM, evidence sources constitute incompatible events, meaning they cannot be considered simultaneously because the joint probability distribution is not defined. In particular, a person cannot represent the 8 dimensional joint distribution of events over disease status (present vs. absent), evidence type (positive vs. negative), and information source (medical history vs. laboratory test). Instead, incompatible events are represented in a lower 4 dimensional space using different bases, which must be evaluated sequentially. Bases correspond to different perspectives of a situation. In the medical diagnosis task, the QOEM assumes that medical history and laboratory tests are viewed one at a time from different perspectives. Order effects arise from this process because the linear algebra operations described below are non-commutative. ","category":"page"},{"location":"model_description/#Basis","page":"Model Description","title":"Basis","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The first step in the process of defining a quantum model is to determine which events are compatible and which are incompatible. The QOEM assumes information source is incompatible, but disease status and evidence type are compatible. Consequentially, the basis of the QOEM consists of four states corresponding to the all possible combinations of disease status and evidence type:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"disease present (p) and positive evidence for disease (p)\ndisease present (p) and negative evidence for disease (n)\ndisease absent (a) and positive evidence for disease (p)\ndisease absent (a) and negative evidence for disease (n)","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"where the values in the parentheses correspond to indices. In the notation used below, the first index corresponds to the presence or absence of the disease, and the second index corresponds to the type of evidence for the disease (positive or negative). Information source is represented as different bases within the reduced 4 dimensional space. Each basis for information source is related to the other bases through a rotation factor, which corresponds to the idea of viewing the diagnosis from different perspectives. More formally, the basis is given by the following orthonormal vectors in the standard position:  ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfB = kettextrmB_ppkettextrmB_pnkettextrmB_apkettextrmB_an","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"where each basis vector consists of a 1 with all other elements equal to zero, e.g., ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"kettextrmB_pp = beginbmatrix\n\t1  \n\t0  \n\t0  \n\t0 \nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Combining the basis vectors into a single matrix, we get the identity matrix:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfI_4 = beginbmatrix\t\t\n\t1  0  0  0\n\t0  1  0  0\n\t0  0  1  0\n\t0  0  0  1\nendbmatrix","category":"page"},{"location":"model_description/#States","page":"Model Description","title":"States","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The state of the cognitive system is a superposition (i.e. linear combination) over basis states:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketboldsymbolpsi = alpha_pp ketB_pp + alpha_pn kettextrmB_pn+ alpha_apkettextrmB_ap+ alpha_an kettextrmB_an","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"where lVertketboldsymbolpsi rVert = 1. The coefficients can be written as:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"boldsymbolalpha = beginbmatrix\n\talpha_pp  \n\talpha_pn  \n\talpha_ap  \n\talpha_an  \nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The initial state is is based on the mean probability judgment after the first symptoms are described: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketboldsymbolpsi = beginbmatrix\n\tsqrt(frac6762)  \n\tsqrt(frac6762)  \n\tsqrt(frac3242)  \n\tsqrt(frac3242)  \nendbmatrix","category":"page"},{"location":"model_description/#Projectors","page":"Model Description","title":"Projectors","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The QOEM makes repeated use of three projection matrices. The following projector is used to evaluate the probability of disease:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfP_d = kettextrmB_pp bratextrmB_pp + kettextrmB_pn bratextrmB_pn = beginbmatrix\t\t\n\t1  0  0  0\n\t0  1  0  0\n\t0  0  0  0\n\t0  0  0  0\nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Notice it spans the 2D sub-space in which the disease is present. The next projector is used to evaluate the probability of positive evidence:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfP_p = kettextrmB_pp bratextrmB_pp + kettextrmB_ap bratextrmB_ap = beginbmatrix\t\t\n\t1  0  0  0\n\t0  0  0  0\n\t0  0  1  0\n\t0  0  0  0\nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Similarly, the projector spans the 2D sub-space in which positive evidence is discovered. Finally, the following projector evaluates the probability of negative evidence:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfP_n = kettextrmB_pn bratextrmB_pn + kettextrmB_an bratextrmB_an = beginbmatrix\t\t\n\t0  0  0  0\n\t0  1  0  0\n\t0  0  0  0\n\t0  0  0  1\nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"As before, the projector spans the 2D sub-space in which negative evidence is discovered. ","category":"page"},{"location":"model_description/#Hamiltonian-Matrices","page":"Model Description","title":"Hamiltonian Matrices","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Hamiltonian matrices govern the decision dynamics of the model. The Hamiltonian matrix mathbfH consists of two components: mathbfH_A is sensitive to the payoff matrix, and mathbfH_B is sensitive to cognitive dissonance between beliefs and actions. The component mathbfH_A is defined as follows: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfH_A = beginbmatrix\t\t\n\tmathbfH_A_d  mathbf0\n\tmathbf0  mathbfH_A_c\nendbmatrix","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"where","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfH_A_k = frac1sqrt1 + mu_k^2beginbmatrix\t\t\n\tmu_k  1\n\t1  -mu_k\nendbmatrix","category":"page"},{"location":"model_description/#Evidence-Evaluation","page":"Model Description","title":"Evidence Evaluation","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"This selection describes the process of selecting an action and determining the defection probability. The time evolution is governed by the unitary transformation matrix which is given by:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfU_h= e^-i cdot t cdot  mathbfH_h","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"mathbfU_l = e^-i cdot t cdot  mathbfH_l","category":"page"},{"location":"model_description/#QOEM-Predictions","page":"Model Description","title":"QOEM Predictions","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"In this section, we go through the steps for computing the predictions in the condition in which the medical history is provided followed by the laboratory test. The predictions follow a similar procedure in the condition in which the order of evidence is reversed.","category":"page"},{"location":"model_description/#Assessment-Without-Evidence-Sources","page":"Model Description","title":"Assessment Without Evidence Sources","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"First, we will compute the probability that the disease is present before additional evidence is collected from medical history and laboratory tests. The computation involves projecting the intial state onto the 2D sub-space spanning disease present and squaring the magnitude of the projection, which is computed as:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Pr(D=d) = lVert mathbfP_d ketpsi rVert^2","category":"page"},{"location":"model_description/#Assessment-With-Medical-History","page":"Model Description","title":"Assessment With Medical History","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Next, suppose we update the assessment of the disease given positive evidence from test i is found. The process involves three computations: (1) projecting the initial state onto the positive evidence sub-space, (2) normalizing the state vector, and (3) projecting the new state onto the 2D sub-space spanning disease present. To update the state, we first project onto the sub-space representing medical history.","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketpsi_p^prime = mathbfP_p mathbfU_hketpsi","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Note that he unitary matrix mathbfU_p changes the system to the basis for medical history. After projecting on to the 2D sub-space for positive, the state collapses onto that sub-space. As a result, the new state must be normalized such that the squared magnitude is 1:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketpsi_p = fracketpsi_p^primelVert ketpsi_p^prime rVert","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The last step involves projecting onto the 2D sub-space spanning disease present and squaring the magnitude to obtain the revised probability judgment for the disease:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Pr(D=d mid S_i = 1) = lVert mathbfP_d ketpsi_p rVert^2","category":"page"},{"location":"model_description/#Assessment-With-Medical-History-and-Laboratory-Test","page":"Model Description","title":"Assessment With Medical History and Laboratory Test","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Next, negative evidence is presented from the laboratory test. A similar sequence of steps is followed. The first step involves projecting onto the 2D sub-space spanned by negative evidence: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketpsi_pn^prime = mathbfP_n mathbfU_l mathbfU_h^dagger ketpsi_p","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Notice that the process of changing bases involves an extra step. Because the model state is currently in the positive medical history basis, it is necessary to perform the inverse operation with the conjugate transpose mathbfU_h^dagger and then switch to the basis for the laboratory test with mathbfU_l. As before, the state collapses to negative evidence because negative evidence was observed. Thus, the state must be normalized:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ketpsi_pn = fracketpsi_pn^primelVert ketpsi_pn^prime rVert","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Now that the state is in the sub-space for negative evidence, the last step involves projecting onto the 2D sub-space for disease present to obtain a probability judgment for the presense of the disease:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Pr(D=d mid S_i = 1 S_j=-1) = lVert mathbfP_d ketpsi_pn rVert^2","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Condition Formula Data Model\n1 Pr(R_1=d mid R_2=d) .84 .81\n2 Pr(R_1=d mid R_2=c) .66 .65\n3 Pr(R_1=d) .55 .57","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The code used to generate the predictions can be viewed by expanding the code block below:","category":"page"},{"location":"model_description/#Dynamics","page":"Model Description","title":"Dynamics","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The plot below shows the dynamics of the model for each condition.","category":"page"},{"location":"model_description/#Interference-Effects","page":"Model Description","title":"Interference Effects","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The plot below shows the interference effect as a function of mu for multiple values of gamma. In the simulations below, we fix t=fracpi2.","category":"page"},{"location":"model_description/#References","page":"Model Description","title":"References","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.","category":"page"},{"location":"basic_usage/#Overview","page":"Basic Usage","title":"Overview","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"This page provides an overview of the API along with examples. ","category":"page"},{"location":"basic_usage/#Make-Predictions","page":"Basic Usage","title":"Make Predictions","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"The quantum order effect model (QOEM) generates predictions for six conditions:","category":"page"},{"location":"basic_usage/#Order-1","page":"Basic Usage","title":"Order 1","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"Pr(disease)\nPr(disease | lab test)\nPr(disease | lab test, medical history)","category":"page"},{"location":"basic_usage/#Order-2","page":"Basic Usage","title":"Order 2","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"Pr(disease)\nPr(disease | medical history)\nPr(disease | medical history, lab test)","category":"page"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"using QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\npredict(model)","category":"page"},{"location":"basic_usage/#Simulate-Model","page":"Basic Usage","title":"Simulate Model","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"The code block below demonstrates how to generate simulated data from the model using rand. In the example, we will generate 100 simulated trials for each condition. ","category":"page"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"using QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nn_trials = 100\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\ndata = rand(model, n_trials)","category":"page"},{"location":"basic_usage/#Evaluate-Log-Likelihood","page":"Basic Usage","title":"Evaluate Log Likelihood","text":"","category":"section"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"The log likelihood of data can be evaluated using logpdf. In the code block below, we generate simulated data and evaluate the logpdf: ","category":"page"},{"location":"basic_usage/","page":"Basic Usage","title":"Basic Usage","text":"using QuantumOrderEffectModels\nΨ = @. √([.35,.35,.15,.15])\nγₕ = 2\nγₗ = .5\nσ = .05\nn_trials = 100\nmodel = QOEM(;Ψ, γₕ, γₗ, σ)\ndata = rand(model, n_trials)\nlogpdf(model, data)","category":"page"},{"location":"#QuantumOrderEffectModels.jl","page":"Home","title":"QuantumOrderEffectModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package contains code for a quantum cognition model of order effects in medical diagonsis.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are two methods for installing the package. Option 1 is to install without version control. In the REPL, use ] to switch to the package mode and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add https://github.com/itsdfish/QuantumOrderEffectModels.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"Option 2 is to install via a custom registry. The advantage of this approach is that you have more control over version control, expecially if you are using a project-specfic environment. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Install the registry using the directions found here.\nAdd the package by typing ] into the REPL and then typing (or pasting):","category":"page"},{"location":"","page":"Home","title":"Home","text":"add QuantumOrderEffectModels","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.","category":"page"}]
}
