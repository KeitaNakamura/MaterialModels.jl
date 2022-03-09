@testset "DruckerPrager" begin
    for elastic in (LinearElastic(E = 1e6, ν = 0.3),)
        for mc_type in (:circumscribed, :inscribed, :plane_strain)
            Random.seed!(1234)
            stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
            yield_function(m, σ) = @matcalc(:yield_function, m; σ)

            ## without tension-cutoff
            m = @inferred DruckerPrager(elastic, mc_type, c = 20.0, ϕ = deg2rad(30), ψ = deg2rad(10), tension_cutoff = false)
            for i in 1:50 # check stress integration by random strain 50 times
                σ = -50.0*rand(SymmetricSecondOrderTensor{3})
                dϵ = @matcalc(:strain, m.elastic; σ)
                σ′ = @inferred stress(m, zero(σ), dϵ)
                f = @inferred yield_function(m, σ)
                if f > 0
                    @test abs(@matcalc(:yield_function, m; σ = σ′)) < sqrt(eps(Float64))
                else
                    @test σ′ ≈ σ
                end
            end

            ## with tension-cutoff
            function recompute_stress(m, σ)
                dϵ = @matcalc(:strain, m.elastic; σ)
                @inferred stress(m, zero(σ), dϵ)
            end
            m′ = @inferred DruckerPrager(elastic, mc_type, c = 20.0, ϕ = deg2rad(30), ψ = deg2rad(10), tension_cutoff = 20.0)

            # zone3 && closed to cap
            σ = 20I + rand(SymmetricSecondOrderTensor{3})
            @test @matcalc(:yield_function, m′; σ) < 0
            @test mean(recompute_stress(m′, σ)) ≈ 20.0

            # zone3 && far from cap
            σ = 50I + rand(SymmetricSecondOrderTensor{3})
            @test @matcalc(:yield_function, m′; σ) > 0
            @test mean(recompute_stress(m′, σ)) ≈ 20.0

            # zone2
            σ = 50*(I + dev(rand(SymmetricSecondOrderTensor{3})))
            @test @matcalc(:yield_function, m′; σ) > 0
            @test mean(recompute_stress(m′, σ)) ≈ 20.0
            @test mean(recompute_stress(m, σ)) > 20.0
            @test abs(@matcalc(:yield_function, m′; σ = recompute_stress(m′, σ))) < sqrt(eps(Float64))
        end
    end
end
