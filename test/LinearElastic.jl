@testset "LinearElastic" begin
    Random.seed!(1234)
    checkparams(x, y) = x.E ≈ y.E && x.K ≈ y.K && x.G ≈ y.G && x.λ ≈ y.λ && x.ν ≈ y.ν
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
    @test @matcalc(:stress, m; σ, dϵ) ≈ σ + @matcalc(:stress, m; ϵ = dϵ)
    @test @matcalc(:stress, m; σ, dϵ) ≈ σ + @matcalc(:stiffness, m) ⊡ dϵ
    @test @matcalc(:stress, m; ϵ = @matcalc(:strain, m; σ)) ≈ σ
end
