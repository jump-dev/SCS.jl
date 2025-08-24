# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    struct ScaledComplexPSDCone <: MOI.AbstractVectorSet
        side_dimension::Int
    end

Similar to `MOI.Scaled{ComplexPositiveSemidefiniteConeTriangle}` but its
vectorization is the lower triangular part column-wise.
"""
struct ScaledComplexPSDCone <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(::Type{ScaledComplexPSDCone}, dim)
    return ScaledComplexPSDCone(isqrt(dim))
end

Base.copy(x::ScaledComplexPSDCone) = ScaledComplexPSDCone(x.side_dimension)

MOI.side_dimension(x::ScaledComplexPSDCone) = x.side_dimension

function MOI.dimension(x::ScaledComplexPSDCone)
    return x.side_dimension^2
end

struct ScaledComplexPSDConeBridge{T,F} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledComplexPSDCone,
    MOI.Scaled{ComplexPositiveSemidefiniteConeTriangle},
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,ScaledComplexPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledComplexPSDConeBridge{T}},
    ::Type{F},
    ::Type{MOI.Scaled{ComplexPositiveSemidefiniteConeTriangle}},
) where {T,F<:MOI.AbstractVectorFunction}
    return ScaledComplexPSDConeBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledComplexPSDConeBridge},
    set::MOI.Scaled{ComplexPositiveSemidefiniteConeTriangle},
)
    return ScaledComplexPSDCone(MOI.side_dimension(set))
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledComplexPSDConeBridge},
    set::ScaledComplexPSDCone,
)
    return MOI.Scaled(
        ComplexPositiveSemidefiniteConeTriangle(set.side_dimension),
    )
end

function _complex_to_scs(func)
    vals = MOI.Utilities.eachscalar(func)
    d = isqrt(length(vals))
    c = 0
    perm = zeros(Int, length(vals))
    for i in 1:d
        c += 1
        perm[c] = i^2
        for j in i+1:d
            triidx = 2i - 1 + (j - 1)^2
            c += 1
            perm[c] = triidx
            c += 1
            perm[c] = triidx+1
        end
    end
    return vals[perm]
end

function _scs_to_complex(func)
    vals = MOI.Utilities.eachscalar(func)
    d = isqrt(length(vals))
    c = 0
    perm = zeros(Int, length(vals))
    for i in 1:d
        c += 1
        perm[i^2] = c
        for j in i+1:d
            triidx = 2i - 1 + (j - 1)^2
            c += 1
            perm[triidx] = c
            c += 1
            perm[triidx+1] = c
        end
    end
    return vals[perm]
end

# Map ConstraintFunction from MOI -> SCS
function MOI.Bridges.map_function(::Type{<:ScaledComplexPSDConeBridge}, f)
    return _complex_to_scs(f)
end

# Used to map the ConstraintPrimal from SCS -> MOI
function MOI.Bridges.inverse_map_function(
    ::Type{<:ScaledComplexPSDConeBridge},
    f,
)
    return _scs_to_complex(f)
end

# Used to map the ConstraintDual from SCS -> MOI
function MOI.Bridges.adjoint_map_function(
    ::Type{<:ScaledComplexPSDConeBridge},
    f,
)
    return _scs_to_complex(f)
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    ::Type{<:ScaledComplexPSDConeBridge},
    f,
)
    return _complex_to_scs(f)
end
