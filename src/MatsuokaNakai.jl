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

function plastic_corrector(model, ϵᵉ_trial33::SymmetricSecondOrderTensor{3})::Tuple{SymmetricSecondOrderTensor{3, Float64}, Bool}
    tol = sqrt(eps(Float64))

    ϵᵉ_trial, m₁, m₂, m₃ = tospectral(ϵᵉ_trial33)
    ϵᵉ = ϵᵉ_trial
    Δγ = 0.0

    for i in 1:20
        σ = @matcalc(:stress, model.elastic; ϵ = ϵᵉ)
        Cᵉ = @matcalc(:principal_compliance, model.elastic)
        dfdσ, f = gradient(σ -> @matcalc(:yieldfunction, model; σ), σ, :all)
        d²gdσ², dgdσ = gradient(σ -> @matcalc(:plasticflow, model; σ), σ, :all)

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
    (model.α * one(SymmetricSecondOrderTensor{3}), false)
end

"""
    @matcalc(:stressall, model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and related variables as `NamedTuple`.
"""
@matcalc_def function stressall(model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    ϵᵉ = @matcalc(:strain, model.elastic; σ)

    # elastic predictor
    ϵᵉ_trial = ϵᵉ + dϵ
    σ_trial = @matcalc(:stress, model.elastic; ϵ = ϵᵉ_trial)
    if @matcalc(:yieldfunction, model; σ = σ_trial) ≤ 0.0
        σ = σ_trial
        return (; σ, status = (plastic = false, converted = true))
    end

    # plastic corrector
    σ, converted = plastic_corrector(model, ϵᵉ_trial)
    (; σ, status = (plastic = true, converted))
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

Compute the yield function.
"""
@matcalc_def function yieldfunction(model::MatsuokaNakai; σ::SymmetricSecondOrderTensor{3})
    α, β = model.α, model.β
    I₁, I₂, I₃ = stress_invariants(σ - α*I)
    -cbrt(I₁*I₂) + cbrt(β*I₃)
end
@matcalc_def function yieldfunction(model::MatsuokaNakai; σ::Vec{3})
    α, β = model.α, model.β
    σ = σ - α*ones(σ)
    I₁ = σ[1] + σ[2] + σ[3]
    I₂ = σ[1]*σ[2] + σ[2]*σ[3] + σ[1]*σ[3]
    I₃ = σ[1] * σ[2] * σ[3]
    -cbrt(I₁*I₂) + cbrt(β*I₃)
end

"""
    @matcalc(:plasticflow, model::MatsuokaNakai; σ::SymmetricSecondOrderTensor{3})

Compute the plastic flow (the gradient of plastic potential function in stress space, i.e., ``\\partial{g}/\\partial{\\sigma}``).
"""
@matcalc_def function plasticflow(model::MatsuokaNakai; σ::SymmetricSecondOrderTensor{3})
    B = model.B
    gradient(σ) do σ
        Base.@_inline_meta
        I₁, I₂, I₃ = stress_invariants(σ)
        -cbrt(I₁*I₂) + cbrt(B*I₃)
    end
end
@matcalc_def function plasticflow(model::MatsuokaNakai; σ::Vec{3})
    B = model.B
    gradient(σ) do σ
        Base.@_inline_meta
        I₁ = σ[1] + σ[2] + σ[3]
        I₂ = σ[1]*σ[2] + σ[2]*σ[3] + σ[1]*σ[3]
        I₃ = σ[1] * σ[2] * σ[3]
        -cbrt(I₁*I₂) + cbrt(B*I₃)
    end
end
