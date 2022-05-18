@testset "MatsuokaNakai" begin
    for elastic in (LinearElastic(E = 1e6, ν = 0.3),)
        Random.seed!(1234)
        stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
        stressall(m, σ, dϵ) = @matcalc(:stressall, m; σ, dϵ)
        yieldfunction(m, σ) = @matcalc(:yieldfunction, m; σ)

        m = @inferred MatsuokaNakai(elastic, c = 20.0, ϕ = deg2rad(30), ψ = deg2rad(10))
        for i in 1:50 # check stress integration by random strain 50 times
            σ = -50.0*rand(SymmetricSecondOrderTensor{3})
            dϵ = @matcalc(:strain, m.elastic; σ)
            ret = @inferred stressall(m, zero(σ), dϵ)
            f = @inferred yieldfunction(m, σ)
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
