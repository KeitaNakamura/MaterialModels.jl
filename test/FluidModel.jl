@testset "NewtonianFluid" begin
    pressure(m, ρ) = @matcalc(:pressure, m; ρ)
    density(m, p) = @matcalc(:density, m; p)
    stress(m, d, ρ) = @matcalc(:stress, m; d, ρ)
    for eos in (MonaghanWaterEOS(ρ_ref = 1e3, B = 1119e3), MorrisWaterEOS(ρ_ref = 1e3, c = 1.5))
        Random.seed!(1234)
        m = NewtonianFluid(eos, μ = 1.01e-3)
        ρ = 1e3 * rand()
        p = @inferred(pressure(m, ρ))
        d = -rand(SymmetricSecondOrderTensor{3})
        @test @inferred(density(m, p)) ≈ ρ
        @test pressure(m, ρ) ≈ pressure(m.eos, ρ)
        @test density(m, p) ≈ density(m.eos, p)
        @test stress(m, d, ρ) ≈ -p*I + m.λ*tr(d)*I + 2*m.μ*d
    end
end

