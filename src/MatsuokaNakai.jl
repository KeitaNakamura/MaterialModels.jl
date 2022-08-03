struct MatsuokaNakai{Elastic <: ElasticModel} <: ElastoPlasticModel
    elastic::Elastic
    α::Float64
    β::Float64
    B::Float64
    c::Float64
    ϕ::Float64
    ψ::Float64
end

function MatsuokaNakai(elastic::Elastic; c::Real, ϕ::Real, ψ::Real = ϕ, checkparameters::Bool=true) where {Elastic <: ElasticModel}
    if checkparameters
        if !(0 ≤ ϕ ≤ 2π)
            @warn "Perhaps you are using degree for internal friction angle `ϕ`, you should use radian instead. This message can be disabled by setting `checkparameters=false` in `MatsuokaNakai` constructor."
        end
        if !(0 ≤ ψ ≤ 2π)
            @warn "Perhaps you are using degree for internal friction angle `ψ`, you should use radian instead. This message can be disabled by setting `checkparameters=false` in `MatsuokaNakai` constructor."
        end
    end
    α = c * cot(ϕ)
    β = (9 - sin(ϕ)^2) / (1 - sin(ϕ)^2)
    B = (9 - sin(ψ)^2) / (1 - sin(ψ)^2)
    MatsuokaNakai{Elastic}(elastic, α, β, B, c, ϕ, ψ)
end

"""
    @matcalc(:stressall, model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and related variables as `NamedTuple`.
"""
@matcalc_def function stressall(model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    α, β = model.α, model.β

    ϵᵉ = @matcalc(:strain, model.elastic; σ)

    # elastic predictor
    ϵᵉ_trial = ϵᵉ + dϵ
    σ_trial = @matcalc(:stress, model.elastic; ϵ = ϵᵉ_trial)

    σ_tip = α * one(σ_trial)
    I₁, I₂, I₃ = stress_invariants(σ_trial - σ_tip)
    if !(I₁ ≤ 0 && I₂ ≥ 0 && I₃ ≤ 0)
        return (; σ = σ_tip, status = (plastic = true, converged = true, tensioncollapse = true))
    end

    if @matcalc(:yieldfunction, model; σ = σ_trial) ≤ 0.0
        σ = σ_trial
        return (; σ, status = (plastic = false, converged = true, tensioncollapse = false))
    end

    # plastic corrector
    σ, converged = plastic_corrector(model, σ_trial; default = σ_tip)
    (; σ, status = (plastic = true, converged, tensioncollapse = !converged))
end

"""
    @matcalc(:stress, model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress.
"""
@matcalc_def function stress(model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    @matcalc(:stressall, model; σ, dϵ).σ
end

"""
    @matcalc(:yieldfunction, model::MatsuokaNakai; σ::SymmetricSecondOrderTensor{3})
    @matcalc(:yieldfunction, model::MatsuokaNakai; σ::Vec{3})

Compute the yield function.
"""
@matcalc_def function yieldfunction(model::MatsuokaNakai; σ::Union{SymmetricSecondOrderTensor{3}, Vec{3}})
    α, β = model.α, model.β
    I₁, I₂, I₃ = stress_invariants(σ - α*delta(σ))
    -cbrt(I₁*I₂) + cbrt(β*I₃)
end

"""
    @matcalc(:plasticflow, model::MatsuokaNakai; σ::SymmetricSecondOrderTensor{3})
    @matcalc(:plasticflow, model::MatsuokaNakai; σ::Vec{3})

Compute the plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plasticflow(model::MatsuokaNakai; σ::Union{SymmetricSecondOrderTensor{3}, Vec{3}})
    B = model.B
    gradient(σ) do σ
        Base.@_inline_meta
        I₁, I₂, I₃ = stress_invariants(σ)
        -cbrt(I₁*I₂) + cbrt(B*I₃)
    end
end
