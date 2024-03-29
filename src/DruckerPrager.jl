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
f = \\| \\bm{s} \\| - (A - Bp)
```

# Plastic flow
```math
g = \\| \\bm{s} \\| + b p
```
"""
function DruckerPrager(elastic::Elastic; A::Real, B::Real, b::Real) where {Elastic <: ElasticModel}
    DruckerPrager{Elastic}(elastic, A, B, b, NaN, NaN, NaN, NaN)
end

"""
    DruckerPrager(::ElasticModel, mohr_coulomb_type; c, ϕ, ψ = ϕ, tensioncutoff = false)

# Parameters
* `mohr_coulomb_type`: choose from `:compression`, `:tension`, `:average`, `:inscribed` and `:planestrain`
* `c`: cohesion
* `ϕ`: internal friction angle (radian)
* `ψ`: dilatancy angle (radian)
* `tensioncutoff`: set limit of mean stress or `false`.
"""
function DruckerPrager(elastic::Elastic, mc_type; c::Real, ϕ::Real, ψ::Real=ϕ, tensioncutoff=false, checkparameters::Bool=true) where {Elastic <: ElasticModel}
    if checkparameters
        if !(0 ≤ ϕ ≤ 2π)
            @warn "Perhaps you are using degree for internal friction angle `ϕ`, you should use radian instead. This message can be disabled by setting `checkparameters=false` in `DruckerPrager` constructor."
        end
        if !(0 ≤ ψ ≤ 2π)
            @warn "Perhaps you are using degree for internal friction angle `ψ`, you should use radian instead. This message can be disabled by setting `checkparameters=false` in `DruckerPrager` constructor."
        end
    end
    mc_type = Symbol(mc_type)
    if mc_type == :compression
        A = 2*√6*c*cos(ϕ) / (3 - sin(ϕ))
        B = 2*√6*sin(ϕ) / (3 - sin(ϕ))
        b = 2*√6*sin(ψ) / (3 - sin(ψ))
    elseif mc_type == :tension
        A = 2*√6*c*cos(ϕ) / (3 + sin(ϕ))
        B = 2*√6*sin(ϕ) / (3 + sin(ϕ))
        b = 2*√6*sin(ψ) / (3 + sin(ψ))
    elseif mc_type == :average
        m_c = DruckerPrager(elastic, :compression; c, ϕ, ψ, tensioncutoff, checkparameters)
        m_t = DruckerPrager(elastic, :tension;     c, ϕ, ψ, tensioncutoff, checkparameters)
        A = (m_c.A + m_t.A) / 2
        B = (m_c.B + m_t.B) / 2
        b = (m_c.b + m_t.b) / 2
    elseif mc_type == :inscribed
        A = 3*√2*c*cos(ϕ) / sqrt(9 + 3sin(ϕ)^2)
        B = 3*√2*sin(ϕ) / sqrt(9 + 3sin(ϕ)^2)
        b = 3*√2*sin(ψ) / sqrt(9 + 3sin(ψ)^2)
    elseif mc_type == :planestrain
        A = 3*√2*c / sqrt(9 + 12tan(ϕ)^2)
        B = 3*√2*tan(ϕ) / sqrt(9 + 12tan(ϕ)^2)
        b = 3*√2*tan(ψ) / sqrt(9 + 12tan(ψ)^2)
    else
        throw(ArgumentError("Choose Mohr-Coulomb type from `:compression`, `:tension`, `:average`, `:inscribed` and `:planestrain`, got $mc_type"))
    end
    if tensioncutoff === true
        throw(ArgumentError("Set value of mean stress limit to enable `tensioncutoff`"))
    elseif tensioncutoff === false
        tensioncutoff = A / B
    elseif tensioncutoff isa Real
        # ok
    else
        throw(ArgumentError("Invalid value for `tensioncutoff`, got $tensioncutoff"))
    end
    DruckerPrager{Elastic}(elastic, A, B, b, c, ϕ, ψ, tensioncutoff)
end

@matcalc_def function stressall(model::DruckerPrager; D_e::SymmetricFourthOrderTensor{3}, σ_trial::SymmetricSecondOrderTensor{3})
    dfdσ, f_trial = gradient(σ_trial -> @matcalc(:yieldfunction, model; σ = σ_trial), σ_trial, :all)
    p_t = model.tensioncutoff
    if f_trial ≤ 0.0 && mean(σ_trial) ≤ p_t
        σ = σ_trial
        return (; σ, status = (plastic = false, tensioncollapse = false))
    end
    # σ, converged = plastic_corrector(model, σ_trial; default = Inf*one(σ_trial))
    dgdσ = @matcalc(:plasticflow, model; σ = σ_trial)
    dλ = f_trial / (dfdσ ⊡ D_e ⊡ dgdσ)
    dϵ_p = dλ * dgdσ
    σ = σ_trial - D_e ⊡ dϵ_p
    if mean(σ) > p_t # should be tension-cutoff
        # very rough tensioncut-off
        # more proper way is to use multi-surface return mapping algorithm?
        #
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
        if @matcalc(:yieldfunction, model; σ) > 0 # zone2 -> find corner
            A, B = model.A, model.B
            p = mean(σ)
            a = A - B*p
            σ = p_t*I + a*normalize(s)
        end
        return (; σ, status = (plastic = true, tensioncollapse = true))
    end
    (; σ, status = (plastic = true, tensioncollapse = false))
end

"""
    @matcalc(:stressall, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and related variables as `NamedTuple`.
"""
@matcalc_def function stressall(model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    D_e = @matcalc(:stiffness, model.elastic)
    σ_trial = @matcalc(:stress, model.elastic; σ, dϵ)
    @matcalc(:stressall, model; D_e, σ_trial)
end

"""
    @matcalc(:stress, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress.
"""
@matcalc_def function stress(model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    @matcalc(:stressall, model; σ, dϵ).σ
end

"""
    @matcalc(:yieldfunction, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the yield function.
"""
@matcalc_def function yieldfunction(model::DruckerPrager; σ::Union{SymmetricSecondOrderTensor{3}, Vec{3}})
    A, B = model.A, model.B
    p = mean(σ)
    s = dev(σ)
    norm(s) - (A - B*p)
end

"""
    @matcalc(:plasticflow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})

Compute the plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plasticflow(model::DruckerPrager; σ::Union{Vec{3}, SymmetricSecondOrderTensor{3}})
    b = model.b
    s = dev(σ)
    s_norm = norm(s)
    if s_norm < sqrt(eps(eltype(σ)))
        dgdσ = b/3 * delta(σ)
    else
        dgdσ = s/s_norm + b/3 * delta(σ)
    end
    dgdσ
end
