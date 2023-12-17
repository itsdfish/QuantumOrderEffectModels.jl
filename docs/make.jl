using Documenter
using QuantumOrderEffectModels

makedocs(
    warnonly = true,
    sitename="QuantumOrderEffectModels",
    format=Documenter.HTML(
        assets=[
            asset(
                "https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap",
                class=:css,
            ),
        ],
        collapselevel=1,
    ),
    modules=[
        QuantumOrderEffectModels, 
        # Base.get_extension(SequentialSamplingModels, :TuringExt),  
        # Base.get_extension(SequentialSamplingModels, :PlotsExt) 
    ],
    pages=[
        "Home" => "index.md",
        "Basic Usage" => "basic_usage.md",
        "Model Description" => "model_description.md",
        "Parameter Estimation" => "parameter_estimation.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo="github.com/itsdfish/QuantumOrderEffectModels.jl.git",
)