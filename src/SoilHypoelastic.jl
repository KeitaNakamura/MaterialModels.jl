struct SoilHypoelastic <: ElasticModel
    κ::Float64
    ν::Float64
    e_0::Float64
    K_p⁻¹::Float64
    G_p⁻¹::Float64
    D_p⁻¹::SymmetricFourthOrderTensor{3, Float64, 36}
end

"""
    SoilHypoelastic(; κ::Real, ν::Real, e_0::Real)

# Equation
```math
\\dot{\\bm{\\sigma}} = 3K \\mathrm{vol}(\\dot{\\bm{\\epsilon}}) + 2G \\mathrm{dev}(\\dot{\\bm{\\epsilon}})
```
where
```math
K = \\frac{1+e_0}{\\kappa} p, \\quad G = \\frac{3(1-\\nu)}{2(1+\\nu)}
```

# Parameters
* `κ`: elastic compressibility index
* `ν`: Poisson's ratio
* `e_0`: initial void ratio
"""
function SoilHypoelastic(; κ::Real, ν::Real, e_0::Real)
    K = (1 + e_0) / κ
    G = 3K * (1-2ν) / 2(1+ν)
    λ = 3K*ν / (1+ν)
    δ = one(SymmetricSecondOrderTensor{3})
    I = one(SymmetricFourthOrderTensor{3})
    D = λ * δ ⊗ δ + 2G * I
    SoilHypoelastic(κ, ν, e_0, K, G, D)
end

"""
    @matcalc(:stress, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress.
"""
@matcalc_def function stress(model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    σ + @matcalc(:stiffness, model; σ) ⊡ dϵ
end

"""
    @matcalc(:stiffness, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})

Compute the 4th-order stiffness tensor.
"""
@matcalc_def function stiffness(model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})
    model.D_p⁻¹ * abs(mean(σ))
end

"""
    @matcalc(:bulkmodulus, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})

Compute the bulk modulus.

```math
\\dot{p} = K \\dot{\\epsilon}_\\mathrm{v}
```
"""
@matcalc_def function bulkmodulus(model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})
    model.K_p⁻¹ * abs(mean(σ))
end

"""
    @matcalc(:shearmodulus, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})

Compute the shear modulus.

```math
\\dot{q} = 3G \\dot{\\epsilon}_\\mathrm{s}
```
"""
@matcalc_def function shearmodulus(model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})
    model.G_p⁻¹ * abs(mean(σ))
end
