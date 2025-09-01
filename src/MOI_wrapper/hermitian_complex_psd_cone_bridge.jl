# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    struct ComplexPositiveSemidefiniteConeTriangle <: MOI.AbstractVectorSet
        side_dimension::Int
    end

Similar to `HermitianPositiveSemidefiniteConeTriangle` but its
vectorization interleaves real and imaginary parts of the off-diagonal.
"""
struct ComplexPositiveSemidefiniteConeTriangle <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(
    ::Type{ComplexPositiveSemidefiniteConeTriangle},
    dim,
)
    return ComplexPositiveSemidefiniteConeTriangle(isqrt(dim))
end

function MOI.side_dimension(x::ComplexPositiveSemidefiniteConeTriangle)
    return x.side_dimension
end

function MOI.dimension(x::ComplexPositiveSemidefiniteConeTriangle)
    return x.side_dimension^2
end

function MOI.Utilities.set_dot(
    x::AbstractVector{S},
    y::AbstractVector{T},
    set::ComplexPositiveSemidefiniteConeTriangle,
) where {S,T}
    U = MA.promote_operation(MA.add_mul, S, T)
    result = zero(U)
    d = set.side_dimension
    k = 0
    for j in 1:d
        for i in 1:j-1
            k += 1
            result = MA.add_mul!!(result, 2, x[k], y[k])
            k += 1
            result = MA.add_mul!!(result, 2, x[k], y[k])
        end
        k += 1
        result = MA.add_mul!!(result, x[k], y[k])
    end
    return result
end

function MOI.Utilities.dot_coefficients(
    a::AbstractVector,
    set::ComplexPositiveSemidefiniteConeTriangle,
)
    d = set.side_dimension
    b = copy(a)
    k = 0
    for j in 1:d
        for i in 1:j-1
            k += 1
            b[k] /= 2
            k += 1
            b[k] /= 2
        end
        k += 1
    end
    return b
end

function MOI.is_set_dot_scaled(::Type{ComplexPositiveSemidefiniteConeTriangle})
    return true
end

struct HermitianComplexPSDConeBridge{T,F} <:
       MOI.Bridges.Constraint.SetMapBridge{
    T,
    ComplexPositiveSemidefiniteConeTriangle,
    MOI.HermitianPositiveSemidefiniteConeTriangle,
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,ComplexPositiveSemidefiniteConeTriangle}
end

MOI.Bridges.bridging_cost(::Type{<:HermitianComplexPSDConeBridge}) = 0.1

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{HermitianComplexPSDConeBridge{T}},
    ::Type{F},
    ::Type{MOI.HermitianPositiveSemidefiniteConeTriangle},
) where {T,F<:MOI.AbstractVectorFunction}
    return HermitianComplexPSDConeBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:HermitianComplexPSDConeBridge},
    set::MOI.HermitianPositiveSemidefiniteConeTriangle,
)
    return ComplexPositiveSemidefiniteConeTriangle(MOI.side_dimension(set))
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:HermitianComplexPSDConeBridge},
    set::ComplexPositiveSemidefiniteConeTriangle,
)
    return MOI.HermitianPositiveSemidefiniteConeTriangle(set.side_dimension)
end

function _hermitian_to_complex(func)
    vals = MOI.Utilities.eachscalar(func)
    dim = length(vals)
    side = isqrt(dim)
    k_re = 1
    k_im = div(side * (side + 1), 2) + 1
    l = 1
    perm = zeros(Int, dim)
    for i in 1:side
        for j in 1:(i-1)
            perm[l] = k_re
            perm[l+1] = k_im
            k_re += 1
            k_im += 1
            l += 2
        end
        perm[l] = k_re
        k_re += 1
        l += 1
    end
    return vals[perm]
end

function _complex_to_hermitian(func)
    vals = MOI.Utilities.eachscalar(func)
    dim = length(vals)
    side = isqrt(dim)
    k_re = 1
    k_im = div(side * (side + 1), 2) + 1
    l = 1
    perm = zeros(Int, dim)
    for i in 1:side
        for j in 1:(i-1)
            perm[k_re] = l
            perm[k_im] = l + 1
            k_re += 1
            k_im += 1
            l += 2
        end
        perm[k_re] = l
        k_re += 1
        l += 1
    end
    return vals[perm]
end

# Map ConstraintFunction from Hermitian -> Complex
function MOI.Bridges.map_function(::Type{<:HermitianComplexPSDConeBridge}, f)
    return _hermitian_to_complex(f)
end

# Used to map the ConstraintPrimal from Complex -> Hermitian
function MOI.Bridges.inverse_map_function(
    ::Type{<:HermitianComplexPSDConeBridge},
    f,
)
    return _complex_to_hermitian(f)
end

# Used to map the ConstraintDual from Complex -> Hermitian
function MOI.Bridges.adjoint_map_function(
    ::Type{<:HermitianComplexPSDConeBridge},
    f,
)
    return _complex_to_hermitian(f)
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    ::Type{<:HermitianComplexPSDConeBridge},
    f,
)
    return _hermitian_to_complex(f)
end
