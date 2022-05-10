struct MatsuokaNakai{Elastic <: ElasticModel} <: ElastoPlasticModel
    elastic::Elastic
    α::Float64
    β::Float64
    B::Float64
    c::Float64
    ϕ::Float64
    ψ::Float64
    tensioncutoff::Float64
end

function MatsuokaNakai(elastic::Elastic; c::Real, ϕ::Real, ψ::Real = ϕ, tensioncutoff=false, checkparameters::Bool=true) where {Elastic <: ElasticModel}
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
    if tensioncutoff === true
        throw(ArgumentError("Set value of mean stress limit to enable `tensioncutoff`"))
    elseif tensioncutoff === false
        tensioncutoff = c * cot(ϕ)
    elseif tensioncutoff isa Real
        # ok
    else
        throw(ArgumentError("Invalid value for `tensioncutoff`, got $tensioncutoff"))
    end
    MatsuokaNakai{Elastic}(elastic, α, β, B, c, ϕ, ψ, tensioncutoff)
end

function compute_stress_on_yieldsurface(model::MatsuokaNakai, p::Real, n::SymmetricSecondOrderTensor{3})
    α = model.α
    β = model.β

    p = p - α

    # coefficients of cubic equation for norm of deviatoric stress
    a = 2β * det(n)
    b = (3-β) * p * (n ⊡ n)
    c = 0
    d = -2(9 - β) * p^3

    # solve cubic equation
    ξ = (-b^2 + 3a*c) / 3a^2
    η = (2b^3 - 9a*b*c + 27a^2*d) / 27a^3
    ζ = sqrt(Complex(27η^2 + 4ξ^3))
    k = 1 / cbrt(6√3)
    s = -b/3a - k * (real((3√3η + ζ)^(1/3)) + real((3√3η - ζ)^(1/3)))

    (p + α)*I + s*n
end

#=
function plastic_corrector(model, σ_trial::SymmetricSecondOrderTensor{3}, ϵᵉ_trial::SymmetricSecondOrderTensor{3})
    tol = sqrt(eps(Float64))

    σ = σ_trial
    ϵᵉ = ϵᵉ_trial
    Δγ = 0.0
    for i in 1:20
        σ = @matcalc(:stress, model.elastic; ϵ = ϵᵉ)
        Cᵉ = @matcalc(:compliance, model.elastic)
        dfdσ, f = gradient(σ -> @matcalc(:yieldfunction, model; σ), σ, :all)
        d²gdσ², dgdσ = gradient(σ -> @matcalc(:plasticflow, model; σ), σ, :all)

        @show f
        R = ϵᵉ - ϵᵉ_trial + Δγ*dgdσ
        norm(R) < tol && abs(f) < tol && return σ

        Ξ = inv(Cᵉ + Δγ * d²gdσ²)
        dfdσ_Ξ = dfdσ ⊡ Ξ
        dΔγ = (f - dfdσ_Ξ ⊡ R) / (dfdσ_Ξ ⊡ dgdσ)
        dϵᵉ = Cᵉ ⊡ Ξ ⊡ (-R - dΔγ * dgdσ)

        Δγ += dΔγ
        ϵᵉ += dϵᵉ
    end
    error()
end
=#

function plastic_corrector(model, σ_trial::SymmetricSecondOrderTensor{3}, ϵᵉ_trial::SymmetricSecondOrderTensor{3})
    tol = sqrt(eps(Float64))

    F = eigen(ϵᵉ_trial)
    n₁ = F.vectors[:,1]
    n₂ = F.vectors[:,2]
    n₃ = F.vectors[:,3]
    m₁ = n₁ ⊗ n₁
    m₂ = n₂ ⊗ n₂
    m₃ = n₃ ⊗ n₃

    ϵᵉ_trial = F.values

    σ = zero(ϵᵉ_trial) # dummy
    ϵᵉ = ϵᵉ_trial
    Δγ = 0.0
    for i in 1:20
        σ = @matcalc(:principal_stress, model.elastic; ϵ = ϵᵉ)
        Cᵉ = @matcalc(:principal_compliance, model.elastic)
        dfdσ, f = gradient(σ -> @matcalc(:yieldfunction, model; σ), σ, :all)
        d²gdσ², dgdσ = gradient(σ -> @matcalc(:plasticflow, model; σ), σ, :all)

        R = ϵᵉ - ϵᵉ_trial + Δγ*dgdσ
        norm(R) < tol && abs(f) < tol && return σ[1]*m₁ + σ[2]*m₂ + σ[3]*m₃

        Ξ = inv(Cᵉ + Δγ * d²gdσ²)
        dfdσ_Ξ = dfdσ ⋅ Ξ
        dΔγ = (f - dfdσ_Ξ ⋅ R) / (dfdσ_Ξ ⋅ dgdσ)
        dϵᵉ = Cᵉ ⋅ Ξ ⋅ (-R - dΔγ * dgdσ)

        Δγ += dΔγ
        ϵᵉ += dϵᵉ
    end
    error()
end

#=
julia> m = MatsuokaNakai(LinearElastic(; E=3e6, ν=0.3333); c = 0.0, ϕ = deg2rad(30.0));
julia> dϵ = SymmetricSecondOrderTensor{3}([-0.5506046120810548 -0.6407529992669043 -0.8385819181737694; -0.6407529992669043 -0.08127479009709193 -0.0006628843765860148; -0.8385819181737694 -0.0006628843765860148 -0.03503232891782471])
julia> @matcalc(:stressall, m; σ = zero(SymmetricSecondOrderTensor{3}), dϵ)
=#

"""
    @matcalc(:stressall, model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute the stress and related variables as `NamedTuple`.
"""
@matcalc_def function stressall(model::MatsuokaNakai{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    ϵᵉ = @matcalc(:strain, model.elastic; σ)
    ϵᵉ_trial = ϵᵉ + dϵ
    σ_trial = @matcalc(:stress, model.elastic; ϵ = ϵᵉ_trial)

    f_trial = @matcalc(:yieldfunction, model; σ = σ_trial)
    p_t = model.tensioncutoff
    if f_trial ≤ 0.0 && mean(σ_trial) ≤ p_t
        σ = σ_trial
        return (; σ, status = (plastic = false, tensioncollapsed = false))
    end

    σ = plastic_corrector(model, σ_trial, ϵᵉ_trial)

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
        if @matcalc(:yieldfunction, model; σ) > 0 # zone2 -> find corner
            σ = compute_stress_on_yieldsurface(model, p_t, normalize(s))
        end
        return (; σ, status = (plastic = true, tensioncollapsed = true))
    end
    (; σ, status = (plastic = true, tensioncollapsed = false))
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
