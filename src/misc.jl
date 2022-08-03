"""
    @matcalc(:soundspeed; K::Real, G::Real, ρ::Real)

Compute the speed of sound in solids.

# Equation
```math
c = \\sqrt{\\frac{K + \\frac{4}{3}G}{\\rho}}
```

# Parameters
* `K`: bulk modulus
* `G`: shear modulus
* `ρ`: density
"""
@matcalc_def function soundspeed(; K::Real, G::Real, ρ::Real)
    T = promote_type(typeof(K), typeof(G), typeof(ρ))
    sqrt((K + T(4/3)*G) / ρ)
end

"""
    @matcalc(:soundspeed; K::Real, G::Real, ρ::Real)

Compute the speed of sound in solids.

# Equation
```math
c = \\sqrt{\\frac{E(1-\\nu)}{\\rho (1+\\nu)(1-2\\nu)}}
```

# Parameters
* `E`: Young's modulus
* `ν`: Poisson's ratio
* `ρ`: density
"""
@matcalc_def function soundspeed(; E::Real, ν::Real, ρ::Real)
    T = promote_type(typeof(E), typeof(ν), typeof(ρ))
    sqrt(E*(1-ν) / (ρ*(1+ν)*(1-2ν)))
end

"""
    @matcalc(:jaumann2caucy; dσ_jaumann::SymmetricSecondOrderTensor{3}, σ::SymmetricSecondOrderTensor{3}, W::SecondOrderTensor{3})

Convert the Jaumann stress rate to the Caucy stress rate.
For example, this function can be used as follows:

```julia
dσᴶ = @matcalc(:stress, model; σ, dϵ = symmetric(∇v*dt)) - σ
dσ = @matcalc(:jaumann2caucy; dσ_jaumann = dσᴶ, σ, W = skew(∇v*dt))
```

# Equation
```math
\\dot{\\bm{\\sigma}}^\\mathrm{J} = \\dot{\\bm{\\sigma}} - \\bm{W} \\cdot \\bm{\\sigma} + \\bm{\\sigma} \\cdot \\bm{W}
```
"""
@matcalc_def function jaumann2caucy(; dσ_jaumann::SymmetricSecondOrderTensor{3}, σ::SymmetricSecondOrderTensor{3}, W::SecondOrderTensor{3})
    dσ_jaumann + symmetric(W ⋅ σ - σ ⋅ W)
end

function tospectral(σ::SymmetricSecondOrderTensor{3})
    F = eigen(σ)
    n₁ = F.vectors[:,1]
    n₂ = F.vectors[:,2]
    n₃ = F.vectors[:,3]
    m₁ = symmetric(n₁ ⊗ n₁, :U)
    m₂ = symmetric(n₂ ⊗ n₂, :U)
    m₃ = symmetric(n₃ ⊗ n₃, :U)
    F.values, m₁, m₂, m₃
end

function fromspectral(σ::Vec{3}, m₁::SymmetricSecondOrderTensor{3}, m₂::SymmetricSecondOrderTensor{3}, m₃::SymmetricSecondOrderTensor{3})
    σ[1]*m₁ + σ[2]*m₂ + σ[3]*m₃
end

function return_mapping(m::MaterialModel, σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3}; default)
    σ_trial, plastic = elastic_predictor(m, σ, dϵ)
    if plastic
        σ, converged = plastic_corrector(m, σ; default)
        (; σ, converged)
    else
        (; σ = σ_trial, converged = true)
    end
end

function elastic_predictor(m::MaterialModel, σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})::Tuple{SymmetricSecondOrderTensor{3, Float64}, Bool}
    σ_trial = @matcalc(:stress, m.elastic; σ, dϵ)
    plastic = @matcalc(:yieldfunction, m; σ = σ_trial) > 0.0
    σ_trial, plastic
end

function plastic_corrector(m::MaterialModel, σ_trial::SymmetricSecondOrderTensor{3}; default)::Tuple{SymmetricSecondOrderTensor{3, Float64}, Bool}
    tol = sqrt(eps(Float64))

    ϵᵉ_trial, m₁, m₂, m₃ = tospectral(@matcalc(:strain, m.elastic; σ = σ_trial))
    ϵᵉ = ϵᵉ_trial
    Δγ = 0.0

    for i in 1:20
        σ = @matcalc(:stress, m.elastic; ϵ = ϵᵉ)
        Cᵉ = @matcalc(:principal_compliance, m.elastic)
        dfdσ, f = gradient(σ -> @matcalc(:yieldfunction, m; σ), σ, :all)
        d²gdσ², dgdσ = gradient(σ -> @matcalc(:plasticflow, m; σ), σ, :all)

        R = ϵᵉ - ϵᵉ_trial + Δγ*dgdσ
        norm(R) < tol && abs(f) < tol && return (fromspectral(σ, m₁, m₂, m₃), true)

        Ξ = inv(Cᵉ + Δγ * d²gdσ²)
        dfdσ_Ξ = dfdσ ⋅ Ξ
        dΔγ = (f - dfdσ_Ξ ⋅ R) / (dfdσ_Ξ ⋅ dgdσ)
        dϵᵉ = Cᵉ ⋅ Ξ ⋅ (-R - dΔγ * dgdσ)

        Δγ += dΔγ
        ϵᵉ += dϵᵉ
    end

    # return tip
    (default, false)
end

# identity unit tensors
delta(x::Type{<: Vec{3}}) = ones(x)
delta(x::Type{<: @Tensor{3,3}}) = one(x)
delta(x::Type{<: @Tensor{@Symmetry{3,3}}}) = one(x)
delta(x::AbstractTensor) = delta(typeof(x))
