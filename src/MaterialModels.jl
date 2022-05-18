module MaterialModels

using Reexport
using MacroTools
@reexport using Tensorial

using Base: @kwdef

export
    @matcalc,
    MaterialModel,
# ElasticModel
    ElasticModel,
    LinearElastic,
    SoilHypoelastic,
# ElastoPlasticModel
    ElastoPlasticModel,
    VonMises,
    DruckerPrager,
    MatsuokaNakai,
# WaterEOS
    WaterEOS,
    MonaghanWaterEOS,
    MorrisWaterEOS,
# FluidModel
    FluidModel,
    NewtonianFluid

abstract type MaterialModel end

include("matcalc.jl")

# elastic models
abstract type ElasticModel <: MaterialModel end
include("LinearElastic.jl")
include("SoilHypoelastic.jl")

# elasto-plastic models
abstract type ElastoPlasticModel <: MaterialModel end
include("VonMises.jl")
include("DruckerPrager.jl")
include("MatsuokaNakai.jl")

# fluids
include("WaterEOS.jl")
include("FluidModel.jl")

include("misc.jl")

end # module
