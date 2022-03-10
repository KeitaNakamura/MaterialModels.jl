abstract type WaterEOS <: MaterialModel end


# https://www.researchgate.net/publication/220789258_Weakly_Compressible_SPH_for_Free_Surface_Flows
"""
    MonaghanWaterEOS(; ρ_ref, B, γ = 7)

# Equation
```math
p = B \\left( \\left( \\frac{\\rho}{\\rho_\\mathrm{ref}} \\right)^\\gamma -1 \\right)
```

# Parameters
* `ρ_ref`: reference density
"""
@kwdef struct MonaghanWaterEOS <: WaterEOS
    ρ_ref::Float64
    B::Float64
    γ::Float64 = 7
end

"""
    @matcalc(:pressure, model::MonaghanWaterEOS; ρ::Real)

Compute the pressure.
"""
@matcalc_def function pressure(model::MonaghanWaterEOS; ρ::Real)
    ρ_ref, B, γ = model.ρ_ref, model.B, model.γ
    B * ((ρ/ρ_ref)^γ - 1)
end

"""
    @matcalc(:density, model::MonaghanWaterEOS; p::Real)

Compute the mass density.
"""
@matcalc_def function density(model::MonaghanWaterEOS; p::Real)
    ρ_ref, B, γ = model.ρ_ref, model.B, model.γ
    ρ_ref * (p/B + 1)^(1/γ)
end


"""
    MorrisWaterEOS(; ρ_ref, c)

# Equation
```math
p = c^2 (\\rho - \\rho_\\mathrm{ref})
```

# Parameters
* `ρ_ref`: reference density
* `c`: speed of sound
"""
@kwdef struct MorrisWaterEOS <: WaterEOS
    ρ_ref::Float64
    c::Float64 # speed of sound
end

"""
    @matcalc(:pressure, model::MorrisWaterEOS; ρ::Real)

Compute the pressure.
"""
@matcalc_def function pressure(model::MorrisWaterEOS; ρ::Real)
    ρ_ref, c = model.ρ_ref, model.c
    c^2 * (ρ - ρ_ref)
end

"""
    @matcalc(:density, model::MorrisWaterEOS; p::Real)

Compute the density.
"""
@matcalc_def function density(model::MorrisWaterEOS; p::Real)
    ρ_ref, c = model.ρ_ref, model.c
    ρ_ref + p/c^2
end
