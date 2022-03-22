"""
    @matcalc(:soundspeed; K::Real, G::Real, ρ::Real)

Compute the speed of sound in solids.

# Equation
```math
c = \\sqrt{\\frac{K + \\frac{4}{3}G}{\\rho}}
```

# Parameters
* `K`: bulk modulus
* `G`: shear modulus
* `ρ`: density
"""
@matcalc_def function soundspeed(; K::Real, G::Real, ρ::Real)
    T = promote_type(typeof(K), typeof(G), typeof(ρ))
    sqrt((K + T(4/3)*G) / ρ)
end

"""
    @matcalc(:soundspeed; K::Real, G::Real, ρ::Real)

Compute the speed of sound in solids.

# Equation
```math
c = \\sqrt{\\frac{E(1-\\nu)}{\\rho (1+\\nu)(1-2\\nu)}}
```

# Parameters
* `E`: Young's modulus
* `ν`: Poisson's ratio
* `ρ`: density
"""
@matcalc_def function soundspeed(; E::Real, ν::Real, ρ::Real)
    T = promote_type(typeof(E), typeof(ν), typeof(ρ))
    sqrt(E*(1-ν) / (ρ*(1+ν)*(1-2ν)))
end

"""
    @matcalc(:jaumann2caucy; dσ_jaumann::SymmetricSecondOrderTensor{3}, σ::SymmetricSecondOrderTensor{3}, W::SecondOrderTensor{3})

Convert the Jaumann stress rate to the Caucy stress rate.
For example, this function can be used as follows:

```julia
dσᴶ = @matcalc(:stress, model; σ, dϵ = symmetric(∇v*dt)) - σ
dσ = @matcalc(:jaumann2caucy; dσ_jaumann = dσᴶ, σ, W = skew(∇v*dt))
```

# Equation
```math
\\dot{\\bm{\\sigma}}^\\mathrm{J} = \\dot{\\bm{\\sigma}} - \\bm{W} \\cdot \\bm{\\sigma} + \\bm{\\sigma} \\cdot \\bm{W}
```
"""
@matcalc_def function jaumann2caucy(; dσ_jaumann::SymmetricSecondOrderTensor{3}, σ::SymmetricSecondOrderTensor{3}, W::SecondOrderTensor{3})
    dσ_jaumann + symmetric(W ⋅ σ - σ ⋅ W)
end
