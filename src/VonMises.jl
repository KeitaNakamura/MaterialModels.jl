struct VonMises{Elastic <: ElasticModel} <: ElastoPlasticModel
    elastic::Elastic
    q_y::Float64
    c::Float64
end

"""
    VonMises(::ElasticModel; q_y)

# Yield function
```math
f = q - q_\\mathrm{y}
```

# Plastic flow
```math
g = q - q_\\mathrm{y}
```
"""
function VonMises(elastic::Elastic; q_y::Real) where {Elastic}
    VonMises{Elastic}(elastic, q_y, NaN)
end

"""
    VonMises(::ElasticModel, mohr_coulomb_type; c)

# Parameters
* `mohr_coulomb_type`: `:plane_strain`
* `c`: cohesion
"""
function VonMises(elastic::Elastic, mc_type; c::Real) where {Elastic}
    mc_type = Symbol(mc_type)
    if mc_type == :plane_strain
        q_y = √3c
    else
        throw(ArgumentError("Only :plane_strain is supported for Mohr-Coulomb type, got $mc_type"))
    end
    VonMises{Elastic}(elastic, q_y, c)
end

"""
    @matcalc(:stress_status, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and also return status `(; plastic::Bool)`.
"""
@matcalc_def function stress_status(model::VonMises{LinearElastic}; D_e::SymmetricFourthOrderTensor{3}, σ_trial::SymmetricSecondOrderTensor{3})
    dfdσ, f_trial = gradient(σ_trial -> @matcalc(:yield_function, model; σ = σ_trial), σ_trial, :all)
    f_trial ≤ 0.0 && return (σ_trial, (plastic = false,))
    dgdσ = @matcalc(:plastic_flow, model; σ = σ_trial)
    Δγ = f_trial / (dgdσ ⊡ D_e ⊡ dfdσ)
    σ = σ_trial - Δγ * (D_e ⊡ dgdσ)
    (σ, (plastic = true,))
end

"""
    @matcalc(:stress_status, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and also return status `(; plastic::Bool)`.
"""
@matcalc_def function stress_status(model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    D_e = @matcalc(:stiffness, model.elastic)
    σ_trial = @matcalc(:stress, model.elastic; σ, dϵ)
    @matcalc(:stress_status, model; D_e, σ_trial)
end

"""
    @matcalc(:stress, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress.
"""
@matcalc_def function stress(model::VonMises; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    first(@matcalc(:stress_status, model; σ, dϵ))
end

"""
    @matcalc(:yield_function, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the yield function.
"""
@matcalc_def function yield_function(model::VonMises; σ::SymmetricSecondOrderTensor{3})
    s = dev(σ)
    q = sqrt(3/2 * s ⊡ s)
    q - model.q_y
end

"""
    @matcalc(:plastic_flow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plastic_flow(model::VonMises; σ::SymmetricSecondOrderTensor{3})
    s = dev(σ)
    _s_ = sqrt(s ⊡ s)
    if _s_ < √eps(eltype(σ))
        dgdσ = zero(s)
    else
        dgdσ = sqrt(3/2) * s / _s_
    end
    dgdσ
end
