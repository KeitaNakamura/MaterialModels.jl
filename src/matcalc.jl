################
# @matcalc_def #
################

join_symbols(kwargs, delim) = Symbol(join(kwargs, delim))

macro matcalc_def(def)
    dict = splitdef(def)
    kwargs = sort(dict[:kwargs]; by = x -> splitarg(x)[1]) # sort by arg name

    dict[:name] = Symbol(:matcalc__, dict[:name], :__, join_symbols(map(x -> splitarg(x)[1], kwargs), :__))
    append!(dict[:args], kwargs)
    empty!(dict[:kwargs])

    esc(combinedef(dict))
end

############
# @matcalc #
############

remove_dots(x::Symbol) = x
remove_dots(x::Expr) = Meta.isexpr(x, :.) ? x.args[2].value : x

macro matcalc(parameters::Expr, val::QuoteNode, model = nothing)
    @assert Meta.isexpr(parameters, :parameters)
    kwargs = sort(parameters.args; by = x -> remove_dots(splitarg(x)[1])) # sort by arg name

    f = Symbol(:matcalc__, val.value, :__, join_symbols(map(x -> remove_dots(splitarg(x)[1]), kwargs), :__))
    args = map(kwargs) do kw
        arg_name, arg_type, slurp, default = splitarg(kw)
        default === nothing ? arg_name : default
    end

    if isnothing(model)
        code = :(MaterialModels.$f($(args...)))
    else
        code = :(MaterialModels.$f($model, $(args...)))
    end
    esc(code)
end

macro matcalc(val::QuoteNode, xs...)
    inds = findall(x -> Meta.isexpr(x, :(=)), xs)
    args = xs[setdiff(1:length(xs), inds)]
    kwargs = xs[inds]
    quote
        @matcalc($val, $(args...); $(kwargs...))
    end |> esc
end

##################
# search_matcalc #
##################

typename(x) = Base.typename(x).name
function _search_matcalc_methods(prefix::Symbol, ::Type{Model}) where {Model <: Union{MaterialModel, Nothing}}
    mod = @__MODULE__
    meths = Method[]
    for name in names(mod, all = true)
        f = getfield(mod, name)
        if isa(f, Base.Callable) && startswith(string(name), string(prefix))
            for m in methods(f, mod)
                if Model <: MaterialModel
                    t = m.sig.parameters[2]
                    if t <: Model || typename(t) == typename(Model)
                        push!(meths, m)
                    end
                else
                    push!(meths, m)
                end
            end
        end
    end
    unique(meths)
end

function search_matcalc_methods(valname::Symbol, Model::Type = Nothing)
    _search_matcalc_methods(Symbol(:matcalc__, valname, :__), Model)
end
function search_matcalc_methods(Model::Type = Nothing)
    _search_matcalc_methods(:matcalc__, Model)
end

struct MatCalc
    method::Method
end

function Base.display(m::MatCalc)
    sig = Tuple{m.method.sig.parameters[2:end]...}
    display(Docs.doc(eval(m.method.name), sig))
end

function Base.show(io::IO, mc::MatCalc)
    m = mc.method
    func_str = string(m.name)
    val = Symbol(only(match(r"matcalc__(.+?)__", func_str).captures))
    args_str = split(replace(func_str, match(r"matcalc__.+?__", func_str).match => ""), "__")
    args_types = m.sig.parameters[3:end]
    str = string("@matcalc(:", val, ", ::", m.sig.parameters[2])
    if args_str != [""]
        str = string(str, "; ", join(map((x...) -> join(x, "::"), args_str, args_types), ", "))
    end
    print(io, str, ")")
end

function search_matcalc(args...)
    map(MatCalc, search_matcalc_methods(args...))
end
