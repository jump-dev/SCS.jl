# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    struct ScaledLogDetConeTriangle <: MOI.AbstractVectorSet
        side_dimension::Int
    end

MOI: [t, u, x] in MOI.Scaled{MOI.LogDetConeTriangle}(n)
SCS: [-t, u, perm(x)] in SCS.ScaledLogDetConeTriangle(n)
"""
struct ScaledLogDetConeTriangle <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(::Type{ScaledLogDetConeTriangle}, dim)
    return ScaledLogDetConeTriangle(div(-1 + isqrt(1 + 8 * dim), 2))
end

Base.copy(set::ScaledLogDetConeTriangle) = ScaledLogDetConeTriangle(set.side_dimension)

function MOI.dimension(set::ScaledLogDetConeTriangle)
    return MOI.dimension(MOI.LogDetConeTriangle(set.side_dimension))
end

struct ScaledLogDetConeTriangleBridge{T,F} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledLogDetConeTriangle,
    MOI.Scaled{MOI.LogDetConeTriangle},
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,ScaledLogDetConeTriangle}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledLogDetConeTriangleBridge{T}},
    ::Type{F},
    ::Type{MOI.Scaled{MOI.LogDetConeTriangle}},
) where {T,F<:MOI.AbstractVectorFunction}
    return ScaledLogDetConeTriangleBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledLogDetConeTriangleBridge},
    set::MOI.Scaled{MOI.LogDetConeTriangle},
)
    return ScaledLogDetConeTriangle(set.set.side_dimension)
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledLogDetConeTriangleBridge},
    set::ScaledLogDetConeTriangle,
)
    return MOI.Scaled(MOI.LogDetConeTriangle(set.side_dimension))
end

function _transform_function(
    ::Type{<:ScaledLogDetConeTriangleBridge{T}},
    func,
    moi_to_scs::Bool,
) where {T}
    fs = MOI.Utilities.eachscalar(func)
    upper_to_lower, lower_to_upper =
        _upper_to_lower_triangular_permutation(length(fs) - 2)
    p = moi_to_scs ? lower_to_upper : upper_to_lower
    return MOI.Utilities.operate(vcat, T, -fs[1], fs[2], fs[p.+2])
    return g
end

# Map ConstraintFunction from MOI -> SCS
function MOI.Bridges.map_function(b::Type{<:ScaledLogDetConeTriangleBridge}, f)
    return _transform_function(b, f, true)
end

# Used to map the ConstraintPrimal from SCS -> MOI
function MOI.Bridges.inverse_map_function(
    b::Type{<:ScaledLogDetConeTriangleBridge},
    f,
)
    return _transform_function(b, f, false)
end

# Used to map the ConstraintDual from SCS -> MOI
function MOI.Bridges.adjoint_map_function(
    b::Type{<:ScaledLogDetConeTriangleBridge},
    f,
)
    return _transform_function(b, f, false)
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    b::Type{<:ScaledLogDetConeTriangleBridge},
    f,
)
    return _transform_function(b, f, true)
end
