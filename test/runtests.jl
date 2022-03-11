using Test
using Random
using MaterialModels

@testset "Utilities" begin
    @testset "@matcalc" begin
        Random.seed!(1234)
        m = LinearElastic(; E = rand(), K = rand())
        σ = rand(SymmetricSecondOrderTensor{3})
        dϵ = rand(SymmetricSecondOrderTensor{3})
        @test @matcalc(:stress, m; σ, dϵ) ≈ @matcalc(:stress, m; dϵ, σ)
    end
    @testset "search_matcalc_methods" begin
        meths = @inferred(MaterialModels.search_matcalc_methods())::Vector{Method}
        @test all(str -> startswith(str, "matcalc__"), map(string, meths))
        meths = @inferred(MaterialModels.search_matcalc_methods(:stress))::Vector{Method}
        @test all(str -> startswith(str, "matcalc__stress__"), map(string, meths))
        # specify model
        meths = @inferred(MaterialModels.search_matcalc_methods(LinearElastic))::Vector{Method}
        @test all(str -> startswith(str, "matcalc__"), map(string, meths))
        @test all(m -> m.sig.parameters[2] <: LinearElastic, meths)
        meths = @inferred(MaterialModels.search_matcalc_methods(:stress, LinearElastic))::Vector{Method}
        @test all(str -> startswith(str, "matcalc__stress__"), map(string, meths))
        @test all(m -> m.sig.parameters[2] <: LinearElastic, meths)
    end
end

# elastic models
include("LinearElastic.jl")

# elasto-plastic models
include("VonMises.jl")
include("DruckerPrager.jl")

# fluids
include("WaterEOS.jl")
