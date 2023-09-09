# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    struct ScaledPSDCone <: MOI.AbstractVectorSet
        side_dimension::Int
    end

Similar to `MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}` but it the
vectorization is the lower triangular part column-wise (or the upper triangular
part row-wise).
"""
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

struct ScaledPSDConeBridge{T,F} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledPSDCone,
    MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,ScaledPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledPSDConeBridge{T}},
    ::Type{F},
    ::Type{MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}},
) where {T,F<:MOI.AbstractVectorFunction}
    return ScaledPSDConeBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
)
    return ScaledPSDCone(MOI.side_dimension(set))
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::ScaledPSDCone,
)
    return MOI.Scaled(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension))
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

function _transform_function(func, moi_to_scs::Bool)
    scalars = MOI.Utilities.eachscalar(func)
    d = length(scalars)
    upper_to_lower, lower_to_upper = _upper_to_lower_triangular_permutation(d)
    if moi_to_scs
        return scalars[lower_to_upper]
    else
        return scalars[upper_to_lower]
    end
end

# Map ConstraintFunction from MOI -> SCS
function MOI.Bridges.map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, true)
end

# Used to map the ConstraintPrimal from SCS -> MOI
function MOI.Bridges.inverse_map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, false)
end

# Used to map the ConstraintDual from SCS -> MOI
function MOI.Bridges.adjoint_map_function(::Type{<:ScaledPSDConeBridge}, f)
    return _transform_function(f, false)
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    ::Type{<:ScaledPSDConeBridge},
    f,
)
    return _transform_function(f, true)
end
