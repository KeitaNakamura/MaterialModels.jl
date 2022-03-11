@testset "VonMises" begin
    for elastic in (LinearElastic(E = 1e6, ν = 0.3),)
        Random.seed!(1234)
        stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
        stress_status(m, σ, dϵ) = @matcalc(:stress_status, m; σ, dϵ)
        yield_function(m, σ) = @matcalc(:yield_function, m; σ)

        ## without tension-cutoff
        m = @inferred VonMises(elastic, :plane_strain, c = 20.0)
        m_drucker = @inferred DruckerPrager(elastic, :plane_strain, ϕ = 0, c = 20.0, tensioncutoff = false)
        for i in 1:50 # check stress integration by random strain 50 times
            σ = -50.0*rand(SymmetricSecondOrderTensor{3})
            dϵ = @matcalc(:strain, m.elastic; σ)
            σ′ = @inferred stress(m, zero(σ), dϵ)
            f = @inferred yield_function(m, σ)
            @test σ′ ≈ stress(m_drucker, zero(σ), dϵ)
            if f > 0
                @test abs(@matcalc(:yield_function, m; σ = σ′)) < sqrt(eps(Float64))
            else
                @test σ′ ≈ σ
            end
        end
    end
end
