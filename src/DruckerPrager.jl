struct DruckerPrager{Elastic <: ElasticModel} <: ElastoPlasticModel
    elastic::Elastic
    A::Float64
    B::Float64
    b::Float64
    tension_cutoff::Float64
end

# for Mohr-Coulomb criterion
function DruckerPrager(elastic::Elastic, mc_type; c::Real, ϕ::Real, ψ::Real = ϕ, tension_cutoff::Union{Real, Bool} = 0) where {Elastic <: ElasticModel}
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
    tension_cutoff === true && throw(ArgumentError("Set value of mean stress limit to enable `tension_cutoff`"))
    if tension_cutoff === false
        tension_cutoff = Inf
    end
    DruckerPrager{Elastic}(elastic, A, B, b, tension_cutoff)
end

@matcalc_def function stress(model::DruckerPrager; D_e::SymmetricFourthOrderTensor{3}, σ_trial::SymmetricSecondOrderTensor{3})
    dfdσ, f_trial = gradient(σ_trial -> @matcalc(:yield_function, model; σ = σ_trial), σ_trial, :all)
    p_t = model.tension_cutoff
    (f_trial ≤ 0.0 && mean(σ_trial) ≤ p_t) && return σ_trial
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
    end
    σ
end

@matcalc_def function stress(model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    D_e = @matcalc(:stiffness, model.elastic)
    σ_trial = @matcalc(:stress, model.elastic; σ, dϵ)
    @matcalc(:stress, model; D_e, σ_trial)
end

@matcalc_def function yield_function(model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})
    A = model.A
    B = model.B
    I₁ = tr(σ)
    s = dev(σ)
    J₂ = tr(s ⋅ s) / 2
    √J₂ - (A - B*I₁)
end

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