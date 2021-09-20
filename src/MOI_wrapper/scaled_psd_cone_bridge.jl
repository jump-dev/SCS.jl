struct ScaledPSDCone <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(::Type{ScaledPSDCone}, dim)
    return ScaledPSDCone(div(-1 + isqrt(1 + 8 * dim), 2))
end

Base.copy(x::ScaledPSDCone) = ScaledPSDCone(x.side_dimension)

MOI.side_dimension(x::ScaledPSDCone) = x.side_dimension

function MOI.dimension(x::ScaledPSDCone)
    return div(x.side_dimension * (x.side_dimension + 1), 2)
end

struct ScaledPSDConeBridge{T,F,G} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledPSDCone,
    MOI.PositiveSemidefiniteConeTriangle,
    F,
    G,
}
    constraint::MOI.ConstraintIndex{F,ScaledPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledPSDConeBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T}
    S = MOI.Utilities.scalar_type(G)
    A = MOI.Utilities.promote_operation(*, T, S, T)
    F = MOI.Utilities.promote_operation(vcat, T, A, S)
    return ScaledPSDConeBridge{T,F,G}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::MOI.PositiveSemidefiniteConeTriangle,
)
    return ScaledPSDCone(set.side_dimension)
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::ScaledPSDCone,
)
    return MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)
end

function _upper_to_lower_triangular_permutation(dim::Int)
    side_dimension = MOI.Utilities.side_dimension_for_vectorized_dimension(dim)
    permutation = zeros(Int, dim)
    i = 0
    for row in 1:side_dimension
        start = div(row * (row + 1), 2)
        for col in row:side_dimension
            i += 1
            permutation[i] = start
            start += col
        end
    end
    return sortperm(permutation), permutation
end

function _transform_function(func, scale, moi_to_scs::Bool)
    d = _output_dimension(func)
    # upper_to_lower[i] maps the i'th element of the upper matrix to the linear
    #   index of the lower
    # lower_to_upper[i] maps the i'th element of the lower matrix to the linear
    #   index of the upper
    upper_to_lower, lower_to_upper = _upper_to_lower_triangular_permutation(d)
    # scale_factor[i] is the amount to scale the i'th linear index in the MOI
    # function by.
    scale_factor = fill(scale, d)
    for i in 1:d
        if MOI.Utilities.is_diagonal_vectorized_index(i)
            scale_factor[i] = one(scale)
        end
    end
    if moi_to_scs
        permutation, inverse_permutation = lower_to_upper, upper_to_lower
    else
        permutation, inverse_permutation = upper_to_lower, lower_to_upper
    end
    return _scaled_permutation(
        func,
        scale_factor,
        permutation,
        inverse_permutation,
        moi_to_scs,
    )
end

# Map ConstraintFunction from MOI -> SCS
function MOI.Bridges.map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, √2, true)
end

# Used to map the ConstraintPrimal from SCS -> MOI
function MOI.Bridges.inverse_map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, 1 / √2, false)
end

# Used to map the ConstraintDual from SCS -> MOI
function MOI.Bridges.adjoint_map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, 1 / √2, false)
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    ::Type{<:ScaledPSDConeBridge},
    f,
)
    return _transform_function(f, √2, true)
end
