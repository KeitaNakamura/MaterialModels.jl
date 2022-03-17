@testset "VonMises" begin
    for elastic in (LinearElastic(E = 1e6, ν = 0.3),)
        Random.seed!(1234)
        stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
        stressall(m, σ, dϵ) = @matcalc(:stressall, m; σ, dϵ)
        yieldfunction(m, σ) = @matcalc(:yieldfunction, m; σ)

        ## without tension-cutoff
        m = @inferred VonMises(elastic, :planestrain, c = 20.0)
        m_drucker = @inferred DruckerPrager(elastic, :planestrain, ϕ = 0, c = 20.0, tensioncutoff = false)
        for i in 1:50 # check stress integration by random strain 50 times
            σ = -50.0*rand(SymmetricSecondOrderTensor{3})
            dϵ = @matcalc(:strain, m.elastic; σ)
            ret = @inferred stressall(m, zero(σ), dϵ)
            f = @inferred yieldfunction(m, σ)
            @test ret.σ ≈ stress(m_drucker, zero(σ), dϵ)
            if f > 0
                @test abs(@matcalc(:yieldfunction, m; ret.σ)) < sqrt(eps(Float64))
                @test ret.status.plastic
            else
                @test ret.σ ≈ σ
                @test !ret.status.plastic
            end
        end
    end
end
