using Documenter
using MaterialModels

makedocs(;
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [MaterialModels],
    sitename = "MaterialModels.jl",
    pages=[
        "Home" => "index.md",
        "Material models" => [
            "LinearElastic.md",
            "VonMises.md",
            "DruckerPrager.md",
            "WaterEOS.md",
            "FluidModel.md",
        ],
    ],
    doctest = true, # :fix
)

deploydocs(
    repo = "github.com/KeitaNakamura/MaterialModels.jl.git",
    devbranch = "main",
)
