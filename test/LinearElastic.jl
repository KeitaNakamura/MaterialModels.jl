@testset "LinearElastic" begin
    Random.seed!(1234)
    checkparams(x, y) = x.E ≈ y.E && x.K ≈ y.K && x.G ≈ y.G && x.λ ≈ y.λ && x.ν ≈ y.ν
    stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
    stress(m, ϵ) = @matcalc(:stress, m; ϵ)
    strain(m, σ) = @matcalc(:strain, m; σ)
    stiffness(m) = @matcalc(:stiffness, m)
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
    @test @inferred(stress(m, @inferred(strain(m, σ)))) ≈ σ
end
