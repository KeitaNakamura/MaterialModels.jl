@testset "WaterEOS" begin
    pressure(m, ρ) = @matcalc(:pressure, m; ρ)
    density(m, p) = @matcalc(:density, m; p)
    for m in (MonaghanWaterEOS(ρ_ref = 1e3, B = 1119e3), MorrisWaterEOS(ρ_ref = 1e3, c = 1.5))
        Random.seed!(1234)
        ρ = 1e3 * rand()
        @test @inferred(density(m, @inferred(pressure(m, ρ)))) ≈ ρ
    end
end
