@testset "LinearElastic" begin
    Random.seed!(1234)
    checkparams(x, y) = x.E ≈ y.E && x.K ≈ y.K && x.G ≈ y.G && x.λ ≈ y.λ && x.ν ≈ y.ν
    stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
    stress(m, ϵ) = @matcalc(:stress, m; ϵ)
    strain(m, σ) = @matcalc(:strain, m; σ)
    stiffness(m) = @matcalc(:stiffness, m)
    compliance(m) = @matcalc(:compliance, m)

    m = LinearElastic(; E = rand(), K = rand())
    @test checkparams(m, LinearElastic(; m.E, m.G))
    @test checkparams(m, LinearElastic(; m.E, m.λ))
    @test checkparams(m, LinearElastic(; m.E, m.ν))
    @test checkparams(m, LinearElastic(; m.K, m.G))
    @test checkparams(m, LinearElastic(; m.K, m.λ))
    @test checkparams(m, LinearElastic(; m.K, m.ν))
    @test checkparams(m, LinearElastic(; m.G, m.λ))
    @test checkparams(m, LinearElastic(; m.G, m.ν))
    @test checkparams(m, LinearElastic(; m.λ, m.ν))

    σ = rand(SymmetricSecondOrderTensor{3})
    dϵ = rand(SymmetricSecondOrderTensor{3})
    @test @inferred(stress(m, σ, dϵ)) ≈ σ + @inferred(stress(m, dϵ))
    @test @inferred(stress(m, σ, dϵ)) ≈ σ + @inferred(stiffness(m)) ⊡ dϵ
    @test @inferred(strain(m, σ)) ≈ @inferred(compliance(m)) ⊡ σ
    @test @inferred(stress(m, @inferred(strain(m, σ)))) ≈ σ

    # principal version
    σ3, m₁, m₂, m₃ = MaterialModels.tospectral(σ)
    dϵ3, m₁′, m₂′, m₃′ = MaterialModels.tospectral(dϵ)
    principal_stiffness(m) = @matcalc(:principal_stiffness, m)
    principal_compliance(m) = @matcalc(:principal_compliance, m)

    @test @inferred(stress(m, σ3, dϵ3)) ≈ σ3 + @inferred(principal_stiffness(m)) ⋅ dϵ3
    @test @inferred(strain(m, σ3)) ≈ @inferred(principal_compliance(m)) ⋅ σ3
    @test MaterialModels.fromspectral(@inferred(stress(m, dϵ3)), m₁′, m₂′, m₃′) ≈ stress(m, dϵ)
    @test MaterialModels.fromspectral(@inferred(strain(m,  σ3)),  m₁,  m₂,  m₃) ≈ strain(m, σ)
end
