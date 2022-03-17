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
* `mohr_coulomb_type`: `:planestrain`
* `c`: cohesion
"""
function VonMises(elastic::Elastic, mc_type; c::Real) where {Elastic}
    mc_type = Symbol(mc_type)
    if mc_type == :planestrain
        q_y = √3c
    else
        throw(ArgumentError("Only :planestrain is supported for Mohr-Coulomb type, got $mc_type"))
    end
    VonMises{Elastic}(elastic, q_y, c)
end

@matcalc_def function stressall(model::VonMises{LinearElastic}; D_e::SymmetricFourthOrderTensor{3}, σ_trial::SymmetricSecondOrderTensor{3})
    dfdσ, f_trial = gradient(σ_trial -> @matcalc(:yieldfunction, model; σ = σ_trial), σ_trial, :all)
    if f_trial ≤ 0.0
        σ = σ_trial
        return (; σ, status = (plastic = false,))
    end
    dgdσ = @matcalc(:plasticflow, model; σ = σ_trial)
    dλ = f_trial / (dgdσ ⊡ D_e ⊡ dfdσ)
    dϵ_p = dλ * dgdσ
    σ = σ_trial - D_e ⊡ dϵ_p
    (; σ, status = (plastic = true,))
end

"""
    @matcalc(:stressall, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and related variables as `NamedTuple`.
"""
@matcalc_def function stressall(model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    D_e = @matcalc(:stiffness, model.elastic)
    σ_trial = @matcalc(:stress, model.elastic; σ, dϵ)
    @matcalc(:stressall, model; D_e, σ_trial)
end

"""
    @matcalc(:stress, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress.
"""
@matcalc_def function stress(model::VonMises; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    @matcalc(:stressall, model; σ, dϵ).σ
end

"""
    @matcalc(:yieldfunction, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the yield function.
"""
@matcalc_def function yieldfunction(model::VonMises; σ::SymmetricSecondOrderTensor{3})
    s = dev(σ)
    q = sqrt(3/2 * s ⊡ s)
    q - model.q_y
end

"""
    @matcalc(:plasticflow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plasticflow(model::VonMises; σ::SymmetricSecondOrderTensor{3})
    s = dev(σ)
    _s_ = sqrt(s ⊡ s)
    if _s_ < √eps(eltype(σ))
        dgdσ = zero(s)
    else
        dgdσ = sqrt(3/2) * s / _s_
    end
    dgdσ
end
