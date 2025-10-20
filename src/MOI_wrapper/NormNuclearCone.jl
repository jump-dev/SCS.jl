# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    struct NormNuclearCone <: MOI.AbstractVectorSet
        row_dim::Int
        column_dim::Int
    end

MOI: [t, X] in MOI.NormNuclearCone(m, n)
SCS: [t, X] in SCS.NormNuclearCone(m, n) where m >= n
"""
struct NormNuclearCone <: MOI.AbstractVectorSet
    row_dim::Int
    column_dim::Int
end

function MOI.dimension(set::NormNuclearCone)
    return set.row_dim * set.column_dim + 1
end

struct NormNuclearConeBridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    set::MOI.NormNuclearCone
    constraint::MOI.ConstraintIndex{MOI.VectorAffineFunction{T},NormNuclearCone}
end

_transpose(::MOI.NormNuclearCone, ::Nothing, ::Bool) = nothing

function _transpose(set::MOI.NormNuclearCone, f::AbstractVector, dir::Bool)
    if set.row_dim >= set.column_dim
        return f
    end
    m, n = set.row_dim, set.column_dim
    return [f[1]; vec(reshape(f[2:end], dir ? (m, n) : (n, m))')]
end

function _transpose(
    set::MOI.NormNuclearCone,
    f::MOI.AbstractVectorFunction,
    dir::Bool,
)
    if set.row_dim >= set.column_dim
        return f
    end
    p = _transpose(set, 1:MOI.output_dimension(f), dir)
    return MOI.Utilities.eachscalar(f)[p]
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{NormNuclearConeBridge{T}},
    model::MOI.ModelLike,
    f::MOI.VectorAffineFunction{T},
    s::MOI.NormNuclearCone,
) where {T}
    m, n = max(s.row_dim, s.column_dim), min(s.row_dim, s.column_dim)
    g = _transpose(s, f, true)
    ci = MOI.add_constraint(model, g, NormNuclearCone(m, n))
    return NormNuclearConeBridge{T}(s, ci)
end

function MOI.supports_constraint(
    ::Type{NormNuclearConeBridge{T}},
    ::Type{MOI.VectorAffineFunction{T}},
    ::Type{MOI.NormNuclearCone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:NormNuclearConeBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{NormNuclearConeBridge{T}},
) where {T}
    return Tuple{Type,Type}[(MOI.VectorAffineFunction{T}, NormNuclearCone)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{NormNuclearConeBridge{T}},
    ::Type{MOI.VectorAffineFunction{T}},
    ::Type{MOI.NormNuclearCone},
) where {T}
    return NormNuclearConeBridge{T}
end

function MOI.get(
    ::NormNuclearConeBridge{T},
    ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T},NormNuclearCone},
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::NormNuclearConeBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{T},NormNuclearCone},
) where {T}
    return [bridge.constraint]
end

function MOI.delete(model::MOI.ModelLike, bridge::NormNuclearConeBridge)
    MOI.delete(model, bridge.constraint)
    return
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::NormNuclearConeBridge,
)
    return bridge.set
end

function MOI.supports(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintPrimalStart,MOI.ConstraintDualStart},
    ::Type{NormNuclearConeBridge{T}},
) where {T}
    return MOI.supports(
        model,
        attr,
        MOI.ConstraintIndex{MOI.VectorAffineFunction{T},NormNuclearCone},
    )
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{
        MOI.ConstraintFunction,
        MOI.ConstraintPrimalStart,
        MOI.ConstraintPrimal,
        MOI.ConstraintDual,
        MOI.ConstraintDualStart,
    },
    bridge::NormNuclearConeBridge{T},
) where {T}
    f = MOI.get(model, attr, bridge.constraint)
    return _transpose(bridge.set, f, false)
end

function MOI.set(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintPrimalStart,MOI.ConstraintDualStart},
    bridge::NormNuclearConeBridge{T},
    start::Union{Nothing,AbstractVector{T}},
) where {T}
    g = _transpose(bridge.set, start, true)
    MOI.set(model, attr, bridge.constraint, start)
    return
end
