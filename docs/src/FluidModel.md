# Models for fluids

## Newtonian fluid

```@docs
NewtonianFluid
```

##### Methods

```@docs
MaterialModels.matcalc__pressure__ρ(::NewtonianFluid, ::Real)
MaterialModels.matcalc__density__p(::NewtonianFluid, ::Real)
MaterialModels.matcalc__stress__d__ρ(::NewtonianFluid, ::SymmetricSecondOrderTensor{3}, ::Real)
```
