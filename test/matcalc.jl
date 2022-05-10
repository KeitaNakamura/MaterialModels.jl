@testset "`matcalc` functions" begin
    @testset "@matcalc" begin
        Random.seed!(1234)
        m = LinearElastic(; E = rand(), K = rand())
        σ = rand(SymmetricSecondOrderTensor{3})
        dϵ = rand(SymmetricSecondOrderTensor{3})
        ρ = rand()
        K, G = m.K, m.G
        @test @matcalc(:stress, m; σ, dϵ) ≈ @matcalc(:stress, m; dϵ, σ)
        @test @matcalc(:stress, m; σ=σ, dϵ=dϵ) ≈ @matcalc(:stress, m; dϵ, σ=σ)
        @test @matcalc(:stress, m, σ=σ, dϵ=dϵ) ≈ @matcalc(:stress, m, dϵ=dϵ, σ=σ)
        @test @matcalc(:soundspeed; K, G, ρ) ≈ @matcalc(:soundspeed; m.G, m.K, ρ)
        @test @matcalc(:soundspeed; K, G=G, ρ) ≈ @matcalc(:soundspeed; G=m.G, m.K, ρ)
        @test @matcalc(:soundspeed, K=K, G=G, ρ=ρ) ≈ @matcalc(:soundspeed, ρ=ρ, K=K, G=G)
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
