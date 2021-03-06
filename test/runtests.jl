using Test
using Random
using MaterialModels

include("matcalc.jl")

# elastic models
include("LinearElastic.jl")
include("SoilHypoelastic.jl")

# elasto-plastic models
include("VonMises.jl")
include("DruckerPrager.jl")
include("MatsuokaNakai.jl")

# fluids
include("WaterEOS.jl")
include("FluidModel.jl")

include("misc.jl")
