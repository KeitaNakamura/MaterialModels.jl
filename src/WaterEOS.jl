abstract type WaterEOS <: MaterialModel end


# https://www.researchgate.net/publication/220789258_Weakly_Compressible_SPH_for_Free_Surface_Flows
@kwdef struct MonaghanWaterEOS <: WaterEOS
    ρ_ref::Float64
    B::Float64
    γ::Float64 = 7
end

@matcalc_def function pressure(model::MonaghanWaterEOS; ρ::Real)
    ρ_ref, B, γ = model.ρ_ref, model.B, model.γ
    B * ((ρ/ρ_ref)^γ - 1)
end

@matcalc_def function density(model::MonaghanWaterEOS; p::Real)
    ρ_ref, B, γ = model.ρ_ref, model.B, model.γ
    ρ_ref * (p/B + 1)^(1/γ)
end


@kwdef struct MorrisWaterEOS <: WaterEOS
    ρ_ref::Float64
    c::Float64 # speed of sound
end

@matcalc_def function pressure(model::MorrisWaterEOS; ρ::Real)
    ρ_ref, c = model.ρ_ref, model.c
    c^2 * (ρ - ρ_ref)
end

@matcalc_def function density(model::MorrisWaterEOS; p::Real)
    ρ_ref, c = model.ρ_ref, model.c
    ρ_ref + p/c^2
end
