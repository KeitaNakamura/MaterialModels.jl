module MaterialModels

using Reexport
using MacroTools
@reexport using Tensorial

export
    MaterialModel,
    LinearElastic,
    DruckerPrager,
    @matcalc


abstract type MaterialModel end
abstract type ElasticModel <: MaterialModel end
abstract type ElastoPlasticModel <: MaterialModel end

################
# @matcalc_def #
################

join_kwargs(kwargs, delim) = Symbol(join(kwargs, delim))

macro matcalc_def(def)
    dict = splitdef(def)
    kwargs = sort(dict[:kwargs]; by = x -> splitarg(x)[1]) # sort by arg name

    dict[:name] = Symbol(:matcalc_, dict[:name], :__, join_kwargs(map(x -> splitarg(x)[1], kwargs), :__))
    append!(dict[:args], kwargs)
    empty!(dict[:kwargs])

    esc(combinedef(dict))
end

############
# @matcalc #
############

macro matcalc(parameters::Expr, val::QuoteNode, model)
    @assert Meta.isexpr(parameters, :parameters)
    kwargs = sort(parameters.args; by = x -> splitarg(x)[1]) # sort by arg name

    f = Symbol(:matcalc_, val.value, :__, join_kwargs(map(x -> splitarg(x)[1], kwargs), :__))
    args = map(kwargs) do kw
        arg_name, arg_type, slurp, default = splitarg(kw)
        default === nothing ? arg_name : default
    end

    quote
        MaterialModels.$f($model, $(args...))
    end |> esc
end

macro matcalc(val::QuoteNode, model)
    f = Symbol(:matcalc_, val.value, :__)
    quote
        MaterialModels.$f($model)
    end |> esc
end

##################
# search_matcalc #
##################

typename(x) = Base.typename(x).name
function _search_matcalc(prefix::Symbol, ::Type{Model}) where {Model <: MaterialModel}
    mod = @__MODULE__
    meths = Method[]
    for name in names(mod, all = true)
        f = getfield(mod, name)
        if isa(f, Base.Callable) && startswith(string(name), string(prefix))
            for m in methods(f, mod)
                t = m.sig.parameters[2]
                if t <: Model || typename(t) == typename(Model)
                    push!(meths, m)
                end
            end
        end
    end
    unique(meths)
end

function search_matcalc(valname::Symbol, ::Type{Model} = MaterialModel) where {Model <: MaterialModel}
    _search_matcalc(Symbol(:matcalc_, valname, :__), Model)
end
function search_matcalc(::Type{Model} = MaterialModel) where {Model <: MaterialModel}
    _search_matcalc(:matcalc_, Model)
end
function search_matcalc(model::MaterialModel)
    _search_matcalc(:matcalc_, typeof(model))
end

include("LinearElastic.jl")
include("DruckerPrager.jl")

end # module
