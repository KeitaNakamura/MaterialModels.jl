@testset "Objective stress rate" begin
    @testset "Jaumann rate" begin
        σ = rand(SymmetricSecondOrderTensor{3})
        σ̇ᴶ = rand(SymmetricSecondOrderTensor{3})
        W = skew(rand(SecondOrderTensor{3}))
        jaumann2caucy(dσ_jaumann, σ, W) = @matcalc(:jaumann2caucy; dσ_jaumann, σ, W)
        σ̇ = @inferred(jaumann2caucy(σ̇ᴶ, σ, W))::SymmetricSecondOrderTensor{3}
        @test σ̇ᴶ ≈ σ̇ - W⋅σ + σ⋅W
    end
end

@testset "Misc" begin
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
