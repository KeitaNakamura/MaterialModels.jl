@testset "SoilHypoelastic" begin
    Random.seed!(1234)
    stress(m, σ, dϵ) = @matcalc(:stress, m; σ, dϵ)
    stiffness(m, σ) = @matcalc(:stiffness, m; σ)
    bulk_modulus(m, σ) = @matcalc(:bulkmodulus, m; σ)
    shear_modulus(m, σ) = @matcalc(:shearmodulus, m; σ)
    mean_deviatoric_stresses(σ) = (mean(σ), √(3/2 * dev(σ) ⊡ dev(σ)))
    m = SoilHypoelastic(; κ = 0.01, ν = 0.3, e_0 = 0.83)
    σ = rand(SymmetricSecondOrderTensor{3})
    dϵ = rand(SymmetricSecondOrderTensor{3})
    p, q = mean_deviatoric_stresses(σ)
    dϵ_v = tr(dϵ)
    dϵ_s = √(2/3 * dev(dϵ) ⊡ dev(dϵ))
    σ′ = @inferred(stress(m, σ, dϵ))
    p′, q′ = mean_deviatoric_stresses(σ′)
    @test gradient(dϵ -> stress(m, σ, dϵ), dϵ) ≈ @inferred(stiffness(m, σ))
    @test p′ ≈ p + @inferred(bulk_modulus(m, σ)) * dϵ_v
    @test q′ ≈ q + 3 * @inferred(shear_modulus(m, σ)) * dϵ_s  rtol=0.005
end
