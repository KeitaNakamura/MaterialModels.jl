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
    @testset "speed of sound" begin
        for T in (Float64, Float32)
            # speed of sound
            Random.seed!(1234)
            m = LinearElastic(; E = rand(), K = rand())
            ρ = rand(T)
            K, G = T(m.K), T(m.G)
            E, ν = T(m.E), T(m.ν)
            soundspeed1(K, G, ρ) = @matcalc(:soundspeed; K, G, ρ)
            soundspeed2(E, ν, ρ) = @matcalc(:soundspeed; E, ν, ρ)
            @test (@inferred soundspeed1(K, G, ρ))::T ≈ (@inferred soundspeed2(E, ν, ρ))::T
        end
    end
end

# elastic models
include("LinearElastic.jl")
include("SoilHypoelastic.jl")

# elasto-plastic models
include("VonMises.jl")
include("DruckerPrager.jl")

# fluids
include("WaterEOS.jl")
include("FluidModel.jl")
