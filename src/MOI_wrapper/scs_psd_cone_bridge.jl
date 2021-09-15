struct SCSPSDCone <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(::Type{SCSPSDCone}, dim)
    return SCSPSDCone(div(-1 + isqrt(1 + 8 * dim), 2))
end

Base.copy(x::SCSPSDCone) = SCSPSDCone(x.side_dimension)

MOI.side_dimension(x::SCSPSDCone) = x.side_dimension

MOI.dimension(x::SCSPSDCone) = div(x.side_dimension * (x.side_dimension + 1), 2)

struct SCSPSDConeBridge{T,F} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    SCSPSDCone,
    MOI.PositiveSemidefiniteConeTriangle,
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,SCSPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SCSPSDConeBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T}
    return SCSPSDConeBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:SCSPSDConeBridge},
    set::MOI.PositiveSemidefiniteConeTriangle,
)
    return SCSPSDCone(set.side_dimension)
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:SCSPSDConeBridge},
    set::SCSPSDCone,
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
    return permutation
end
function MOI.Bridges.map_function(
    ::Type{<:SCSPSDConeBridge{T}},
    func,
) where {T}
    scalars = MOI.Utilities.eachscalar(func)
    p = _upper_to_lower_triangular_permutation(length(scalars))
    return scalars[p]
end

function MOI.Bridges.inverse_map_function(B::Type{<:SCSPSDConeBridge}, f)
    return MOI.Bridges.map_function(B, f)
end

function MOI.Bridges.adjoint_map_function(B::Type{<:SCSPSDConeBridge}, f)
    return MOI.Bridges.map_function(B, f)
end

function MOI.Bridges.inverse_adjoint_map_function(
    B::Type{<:SCSPSDConeBridge},
    f,
)
    return MOI.Bridges.map_function(B, f)
end
