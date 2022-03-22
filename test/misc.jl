@testset "Objective stress rate" begin
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
