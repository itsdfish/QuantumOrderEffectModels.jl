module QuantumOrderEffectModels

    using Distributions: Binomial
    using Distributions: ContinuousUnivariateDistribution
    using LinearAlgebra

    import Distributions: logpdf
    import Distributions: pdf
    import Distributions: rand 
    import Distributions: loglikelihood

    export predict
    export loglikelihood
    export logpdf 
    export pdf
    export QOEM
    export rand 

    include("structs.jl")
    include("functions.jl")

end
