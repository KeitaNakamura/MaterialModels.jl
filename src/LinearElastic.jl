struct LinearElastic <: ElasticModel
    E::Float64
    K::Float64
    G::Float64
    λ::Float64
    ν::Float64
    D::SymmetricFourthOrderTensor{3, Float64, 36}
    Dinv::SymmetricFourthOrderTensor{3, Float64, 36}
    D33::Mat{3, 3, Float64, 9}
    D33inv::Mat{3, 3, Float64, 9}
end

"""
    LinearElastic(; parameters...)

# Parameters
Choose only 2 parameters.
* `E`: Young's modulus
* `K`: bulk modulus
* `G`: shear modulus
* `λ`: Lamé's first parameter
* `ν`: Poisson's ratio
"""
function LinearElastic(; kwargs...)
    @assert length(kwargs) == 2
    params = values(kwargs)
    if haskey(params, :K)
        K = params.K
        if haskey(params, :E)
            E = params.E
            λ = 3K*(3K-E) / (9K-E)
            G = 3K*E / (9K-E)
            ν = (3K-E) / 6K
        elseif haskey(params, :λ)
            λ = params.λ
            E = 9K*(K-λ) / (3K-λ)
            G = 3(K-λ) / 2
            ν = λ / (3K-λ)
        elseif haskey(params, :G)
            G = params.G
            E = 9K*G / (3K+G)
            λ = K - 2G/3
            ν = (3K-2G) / 2(3K+G)
        elseif haskey(params, :ν)
            ν = params.ν
            E = 3K*(1-2ν)
            λ = 3K*ν / (1+ν)
            G = 3K*(1-2ν) / 2(1+ν)
        end
    elseif haskey(params, :E)
        E = params.E
        if haskey(params, :λ)
            λ = params.λ
            R = √(E^2 + 9λ^2 + 2E*λ)
            K = (E+3λ+R) / 6
            G = (E-3λ+R) / 4
            ν = 2λ / (E+λ+R)
        elseif haskey(params, :G)
            G = params.G
            K = E*G / 3(3G-E)
            λ = G*(E-2G) / (3G-E)
            ν = E/2G - 1
        elseif haskey(params, :ν)
            ν = params.ν
            K = E / 3(1-2ν)
            λ = E*ν / ((1+ν)*(1-2ν))
            G = E / 2(1+ν)
        end
    elseif haskey(params, :λ)
        λ = params.λ
        if haskey(params, :G)
            G = params.G
            K = λ + 2G/3
            E = G*(3λ+2G) / (λ+G)
            ν = λ / 2(λ+G)
        elseif haskey(params, :ν)
            ν = params.ν
            K = λ*(1+ν) / 3ν
            E = λ*(1+ν)*(1-2ν) / ν
            G = λ*(1-2ν) / 2ν
        end
    elseif haskey(params, :G)
        G = params.G
        if haskey(params, :ν)
            ν = params.ν
            K = 2G*(1+ν) / 3(1-2ν)
            E = 2G*(1+ν)
            λ = 2G*ν / (1-2ν)
        end
    end

    δ = one(SymmetricSecondOrderTensor{3})
    I = one(SymmetricFourthOrderTensor{3})
    D = λ * δ ⊗ δ + 2G * I

    δ = ones(Vec{3})
    I = one(Mat{3,3})
    D33 = λ * δ ⊗ δ + 2G * I

    LinearElastic(E, K, G, λ, ν, D, inv(D), D33, inv(D33))
end

"""
    @matcalc(:stress, model::LinearElastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})

Compute stress.
"""
@matcalc_def function stress(model::LinearElastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})
    σ + model.D ⊡ dϵ
end

"""
    @matcalc(:stress, model::LinearElastic; ϵ::SymmetricSecondOrderTensor{3})

Compute stress.
"""
@matcalc_def function stress(model::LinearElastic; ϵ::SymmetricSecondOrderTensor{3})
    model.D ⊡ ϵ
end

"""
    @matcalc(:strain, model::LinearElastic; σ::SymmetricSecondOrderTensor{3})

Compute strain.
"""
@matcalc_def function strain(model::LinearElastic; σ::SymmetricSecondOrderTensor{3})
    model.Dinv ⊡ σ
end

"""
    @matcalc(:stiffness, model::LinearElastic)

Return fourth-order stiffness tensor.
"""
@matcalc_def function stiffness(model::LinearElastic)
    model.D
end

"""
    @matcalc(:compliance, model::LinearElastic)

Return fourth-order compliance tensor.
"""
@matcalc_def function compliance(model::LinearElastic)
    model.Dinv
end


# principal versions

"""
    @matcalc(:stress, model::LinearElastic; σ::Vec{3}, dϵ::Vec{3})

Compute principal stresses.
"""
@matcalc_def function stress(model::LinearElastic; σ::Vec{3}, dϵ::Vec{3})
    σ + model.D33 ⋅ dϵ
end

"""
    @matcalc(:stress, model::LinearElastic; ϵ::Vec{3})

Compute principal stresses.
"""
@matcalc_def function stress(model::LinearElastic; ϵ::Vec{3})
    model.D33 ⋅ ϵ
end

"""
    @matcalc(:strain, model::LinearElastic; σ::Vec{3})

Compute principal strains.
"""
@matcalc_def function strain(model::LinearElastic; σ::Vec{3})
    model.D33inv ⋅ σ
end

"""
    @matcalc(:principal_stiffness, model::LinearElastic)

Return stiffness tensor in principal stress space.
"""
@matcalc_def function principal_stiffness(model::LinearElastic)
    model.D33
end

"""
    @matcalc(:principal_compliance, model::LinearElastic)

Return compliance tensor in principal stress space.
"""
@matcalc_def function principal_compliance(model::LinearElastic)
    model.D33inv
end
