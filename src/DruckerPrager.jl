struct DruckerPrager{Elastic <: ElasticModel} <: ElastoPlasticModel
    elastic::Elastic
    A::Float64
    B::Float64
    b::Float64
    c::Float64
    ϕ::Float64
    ψ::Float64
    tensioncutoff::Float64
end

"""
    DruckerPrager(::ElasticModel; A, B, b)

# Yield function
```math
f = \\sqrt{J_2} - (A - B I_1)
```

# Plastic flow
```math
g = \\sqrt{J_2} + b I_1
```
"""
function DruckerPrager(elastic::Elastic; A::Real, B::Real, b::Real) where {Elastic <: ElasticModel}
    DruckerPrager{Elastic}(elastic, A, B, b, NaN, NaN, NaN, NaN)
end

"""
    DruckerPrager(::ElasticModel, mohr_coulomb_type; c, ϕ, ψ = ϕ, tensioncutoff = 0)

# Parameters
* `mohr_coulomb_type`: choose from `:circumscribed`, `:inscribed` and `:plane_strain`
* `c`: cohesion
* `ϕ`: internal friction angle
* `ψ`: dilatancy angle
* `tensioncutoff`: set limit of mean stress or `false`
"""
function DruckerPrager(elastic::Elastic, mc_type; c::Real, ϕ::Real, ψ::Real = ϕ, tensioncutoff::Union{Real, Bool} = 0) where {Elastic <: ElasticModel}
    mc_type = Symbol(mc_type)
    if mc_type == :circumscribed
        A = 6c*cos(ϕ) / (√3 * (3 - sin(ϕ)))
        B = 2sin(ϕ) / (√3 * (3 - sin(ϕ)))
        b = 2sin(ψ) / (√3 * (3 - sin(ψ)))
    elseif mc_type == :inscribed
        A = 6c*cos(ϕ) / (√3 * (3 + sin(ϕ)))
        B = 2sin(ϕ) / (√3 * (3 + sin(ϕ)))
        b = 2sin(ψ) / (√3 * (3 + sin(ψ)))
    elseif mc_type == :plane_strain
        A = 3c / sqrt(9 + 12tan(ϕ)^2)
        B = tan(ϕ) / sqrt(9 + 12tan(ϕ)^2)
        b = tan(ψ) / sqrt(9 + 12tan(ψ)^2)
    else
        throw(ArgumentError("Choose Mohr-Coulomb type from :circumscribed, :inscribed and :plane_strain, got $mc_type"))
    end
    tensioncutoff === true && throw(ArgumentError("Set value of mean stress limit to enable `tensioncutoff`"))
    if tensioncutoff === false
        tensioncutoff = Inf
    end
    DruckerPrager{Elastic}(elastic, A, B, b, c, ϕ, ψ, tensioncutoff)
end

@matcalc_def function stress_status(model::DruckerPrager; D_e::SymmetricFourthOrderTensor{3}, σ_trial::SymmetricSecondOrderTensor{3})
    dfdσ, f_trial = gradient(σ_trial -> @matcalc(:yield_function, model; σ = σ_trial), σ_trial, :all)
    p_t = model.tensioncutoff
    (f_trial ≤ 0.0 && mean(σ_trial) ≤ p_t) && return (σ_trial, (plastic = false, tensioncutoff = false))
    dgdσ = @matcalc(:plastic_flow, model; σ = σ_trial)
    Δγ = f_trial / (dgdσ ⊡ D_e ⊡ dfdσ)
    σ = σ_trial - Δγ * (D_e ⊡ dgdσ)
    if mean(σ) > p_t # should be tension-cutoff
        # \<- yield surface
        #  \
        #   \          (1)
        #    \   <-------------*
        #     \                      zone2
        #      \ corner ____________________
        #       |
        # zone1 |      (1)           zone3
        #       |<-------------*
        #       |
        # ---------------> p
        s = dev(σ_trial) # NOT `σ`
        σ = p_t*I + s # slide stress along p axis (process (1))
        if @matcalc(:yield_function, model; σ) > 0 # zone2 -> find corner
            A = model.A
            B = model.B
            I₁ = tr(σ)
            J₂ = (A - B*I₁)^2
            a = sqrt(2J₂ / (s ⊡ s))
            σ = p_t*I + a*s
        end
        return (σ, (plastic = true, tensioncutoff = true))
    end
    (σ, (plastic = true, tensioncutoff = false))
end

"""
    @matcalc(:stress_status, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute stress and also return status `(; plastic::Bool, tensioncutoff::Bool)`.
"""
@matcalc_def function stress_status(model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    D_e = @matcalc(:stiffness, model.elastic)
    σ_trial = @matcalc(:stress, model.elastic; σ, dϵ)
    @matcalc(:stress_status, model; D_e, σ_trial)
end

"""
    @matcalc(:stress, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute stress.
"""
@matcalc_def function stress(model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    first(@matcalc(:stress_status, model; σ, dϵ))
end

"""
    @matcalc(:yield_function, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute yield function.
"""
@matcalc_def function yield_function(model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})
    A = model.A
    B = model.B
    I₁ = tr(σ)
    s = dev(σ)
    J₂ = tr(s ⋅ s) / 2
    √J₂ - (A - B*I₁)
end

"""
    @matcalc(:plastic_flow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plastic_flow(model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})
    b = model.b
    s = dev(σ)
    J₂ = tr(s ⋅ s) / 2
    if J₂ < eps(typeof(J₂))
        dgdσ = b * one(σ)
    else
        dgdσ = s / (2*√J₂) + b * one(σ)
    end
    dgdσ
end
