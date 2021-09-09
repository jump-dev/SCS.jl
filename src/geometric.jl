mutable struct GeometricConicForm{T,AT,C} <: MOI.ModelLike
    cone_types::C
    cone_types_dict::Dict{DataType,Int}
    num_rows::Vector{Int}
    dimension::Dict{Int,Int}
    A::Union{Nothing,AT} # The constraints are
    b::Vector{T}          # `b - Ax in cones`
    sense::MOI.OptimizationSense
    objective_constant::T # The objective is
    c::Vector{T}          # `sense c'x + objective_constant`
    primal::Vector{T}
    slack::Vector{T}
    dual::Vector{T}
    function GeometricConicForm{T,AT}(cone_types) where {T,AT}
        model = new{T,AT,typeof(cone_types)}()
        model.cone_types = cone_types
        model.cone_types_dict =
            Dict{DataType,Int}(s => i for (i, s) in enumerate(cone_types))
        model.num_rows = zeros(Int, length(cone_types))
        model.dimension = Dict{Int,Int}()
        model.A = nothing
        model.primal = T[]
        model.slack = T[]
        model.dual = T[]
        return model
    end
end

MOI.is_empty(model::GeometricConicForm) = model.A === nothing

function MOI.empty!(model::GeometricConicForm{T}) where {T}
    empty!(model.dimension)
    fill!(model.num_rows, 0)
    model.A = nothing
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective_constant = zero(T)
    return
end

function MOI.supports_constraint(
    model::GeometricConicForm,
    ::Type{MOI.VectorAffineFunction{Float64}},
    ::Type{S},
) where {S<:MOI.AbstractVectorSet}
    return haskey(model.cone_types_dict, S)
end

function _allocate_variables(
    model::GeometricConicForm{T,AT},
    vis_src,
    idxmap,
) where {T,AT}
    model.A = AT(length(vis_src))
    for (i, vi) in enumerate(vis_src)
        idxmap[vi] = MOI.VariableIndex(i)
    end
    return
end

function rows(
    model::GeometricConicForm,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
)
    return ci.value .+ (1:model.dimension[ci.value])
end

function MOI.supports(
    ::GeometricConicForm,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.set(
    ::GeometricConicForm,
    ::MOI.VariablePrimalStart,
    ::MOI.VariableIndex,
    ::Nothing,
)
    return
end

function MOI.set(
    model::GeometricConicForm,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Float64,
)
    model.primal[vi.value] = value
    return
end

function MOI.supports(
    ::GeometricConicForm,
    ::MOI.ConstraintPrimalStart,
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

function MOI.set(
    ::GeometricConicForm,
    ::MOI.ConstraintPrimalStart,
    ::MOI.ConstraintIndex,
    ::Nothing,
)
    return
end

function MOI.set(
    model::GeometricConicForm,
    ::MOI.ConstraintPrimalStart,
    ci::MOI.ConstraintIndex,
    value,
)
    # TODO(odow): why is offset not used?
    offset = constroffset(model, ci)
    model.slack[rows(model, ci)] .= value
    return
end

function MOI.supports(
    ::GeometricConicForm,
    ::MOI.ConstraintDualStart,
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

function MOI.set(
    ::GeometricConicForm,
    ::MOI.ConstraintDualStart,
    ::MOI.ConstraintIndex,
    ::Nothing,
)
    return
end

function MOI.set(
    model::GeometricConicForm,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex,
    value,
)
    # TODO(odow): why is offset not used?
    offset = constroffset(model, ci)
    model.dual[rows(model, ci)] .= value
    return
end

function MOI.set(
    model::GeometricConicForm,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.sense = sense
    return
end

variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value

function variable_index_value(t::MOI.VectorAffineTerm)
    return variable_index_value(t.scalar_term)
end

function MOI.set(
    model::GeometricConicForm,
    ::MOI.ObjectiveFunction,
    f::MOI.ScalarAffineFunction{Float64},
)
    c = Vector(
        sparsevec(
            variable_index_value.(f.terms),
            MOI.coefficient.(f.terms),
            model.A.n,
        ),
    )
    model.objective_constant = f.constant
    model.c = c
    return
end

function _allocate_constraint(
    model::GeometricConicForm,
    src,
    indexmap,
    cone_id,
    ci,
)
    # TODO use `CanonicalConstraintFunction`
    func = MOI.get(src, MOI.ConstraintFunction(), ci)
    func =
        MOI.Utilities.is_canonical(func) ? func : MOI.Utilities.canonical(func)
    allocate_terms(model.A, indexmap, func)
    offset = model.num_rows[cone_id]
    model.num_rows[cone_id] = offset + MOI.output_dimension(func)
    return ci, offset, func
end

function _allocate_constraints(
    model::GeometricConicForm,
    src,
    indexmap,
    cone_id,
    ::Type{S},
) where {S}
    cis = MOI.get(
        src,
        MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64},S}(),
    )
    return map(cis) do ci
        return _allocate_constraint(model, src, indexmap, cone_id, ci)
    end
end

function _load_variables(model::GeometricConicForm, nvars::Integer)
    m = sum(model.num_rows)
    model.A.m = m
    model.b = zeros(m)
    model.c = zeros(model.A.n)
    allocate_nonzeros(model.A)
    # `model.primal` contains the result of the previous optimization.
    # It is used as a warm start if its length is the same, e.g.
    # probably because no variable and/or constraint has been added.
    if length(model.primal) != nvars
        model.primal = zeros(nvars)
    end
    @assert length(model.dual) == length(model.slack)
    if length(model.dual) != m
        model.dual = zeros(m)
        model.slack = zeros(m)
    end
    return
end

function _load_constraints(
    model::GeometricConicForm,
    src,
    indexmap,
    cone_offset,
    i,
    cache,
    preprocess,
)
    for (ci_src, offset_in_cone, func) in cache
        offset = cone_offset + offset_in_cone
        set = MOI.get(src, MOI.ConstraintSet(), ci_src)
        new_func = preprocess(func, set)
        load_terms(model.A, indexmap, new_func, offset)
        copyto!(model.b, offset + 1, new_func.constants)
        model.dimension[offset] = MOI.output_dimension(func)
        indexmap[ci_src] = typeof(ci_src)(offset)
    end
end

preprocess(func, set) = func

function MOI.copy_to(
    dest::GeometricConicForm,
    src::MOI.ModelLike;
    preprocess = preprocess,
    copy_names::Bool = true,
)
    MOI.empty!(dest)

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    idxmap = MOI.Utilities.IndexMap()

    has_constraints = BitSet()
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        i = get(dest.cone_types_dict, S, nothing)
        if i === nothing || F != MOI.VectorAffineFunction{Float64}
            throw(MOI.UnsupportedConstraint{F,S}())
        end
        push!(has_constraints, i)
    end

    _allocate_variables(dest, vis_src, idxmap)

    # Allocate constraints
    caches = map(collect(has_constraints)) do i
        return _allocate_constraints(dest, src, idxmap, i, dest.cone_types[i])
    end

    # Load variables
    _load_variables(dest, length(vis_src))

    # Set variable attributes
    MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap, vis_src)

    # Set model attributes
    MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap)

    # Load constraints
    offset = 0
    for (i, cache) in zip(has_constraints, caches)
        _load_constraints(dest, src, idxmap, offset, i, cache, preprocess)
        offset += dest.num_rows[i]
    end

    final_touch(dest.A)

    return idxmap
end
