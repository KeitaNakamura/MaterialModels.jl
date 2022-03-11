abstract type FluidModel <: MaterialModel end

struct NewtonianFluid{EOS <: WaterEOS} <: FluidModel
    eos::EOS
    μ::Float64
    λ::Float64
end

"""
    NewtonianFluid(::WaterEOS; μ, λ = -2μ/3)

# Equation
```math
\\sigma_{ij} = -p \\delta_{ij} + \\lambda d_{kk} \\delta_{ij} + 2\\mu d_{ij},
```

where ``p`` is the pressure and ``d_{ij}`` is the rate of deformation tensor.

# Parameters
* `μ`: dynamic viscosity
* `λ`: second coefficient of viscosity
"""
NewtonianFluid(eos::EOS; μ::Real, λ::Real = -2μ/3) where {EOS} = NewtonianFluid{EOS}(eos, μ, λ)

"""
    @matcalc(:pressure, model::NewtonianFluid{EOS}; ρ::Real)

Compute the pressure based on the equation of state `EOS`.
"""
@matcalc_def function pressure(m::NewtonianFluid; ρ::Real)
    @matcalc(:pressure, m.eos; ρ)
end

"""
    @matcalc(:pressure, model::NewtonianFluid{EOS}; p::Real)

Compute the mass density based on the equation of state `EOS`.
"""
@matcalc_def function density(m::NewtonianFluid; p::Real)
    @matcalc(:density, m.eos; p)
end

"""
    @matcalc(:stress, model::NewtonianFluid; d::SymmetricSecondOrderTensor{3}, ρ::Real)

Compute the stress.
"""
@matcalc_def function stress(m::NewtonianFluid; d::SymmetricSecondOrderTensor{3}, ρ::Real)
    μ, λ = m.μ, m.λ
    p = @matcalc(:pressure, m; ρ)
    -p*I + λ*tr(d)*I + 2μ*d
end
