# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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

struct ScaledPSDConeBridge{T,G} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledPSDCone,
    MOI.PositiveSemidefiniteConeTriangle,
    MOI.VectorAffineFunction{T},
    G,
}
    constraint::MOI.ConstraintIndex{MOI.VectorAffineFunction{T},ScaledPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledPSDConeBridge{T}},
    ::Type{G},
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T,G<:Union{MOI.VectorOfVariables,MOI.VectorAffineFunction{T}}}
    return ScaledPSDConeBridge{T,G}
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

function _transform_function(
    func::MOI.VectorAffineFunction{T},
    scale,
    moi_to_scs::Bool,
) where {T}
    d = MOI.output_dimension(func)
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
            scale_factor[i] = 1.0
        end
    end
    scaled_constants = if moi_to_scs
        (func.constants.*scale_factor)[lower_to_upper]
    else
        func.constants[upper_to_lower] .* scale_factor
    end
    scaled_terms = MOI.VectorAffineTerm{T}[]
    for term in func.terms
        row = term.output_index
        i = moi_to_scs ? upper_to_lower[row] : lower_to_upper[row]
        scale_i = moi_to_scs ? scale_factor[row] : scale_factor[i]
        push!(
            scaled_terms,
            MOI.VectorAffineTerm(
                i,
                MOI.ScalarAffineTerm(
                    term.scalar_term.coefficient * scale_i,
                    term.scalar_term.variable,
                ),
            ),
        )
    end
    return MOI.VectorAffineFunction(scaled_terms, scaled_constants)
end

function _transform_function(func::MOI.VectorOfVariables, scale, moi_to_scs)
    new_f = MOI.Utilities.operate(*, Float64, 1.0, func)
    return _transform_function(new_f, scale, moi_to_scs)
end

function _transform_function(func::Vector{T}, scale, moi_to_scs::Bool) where {T}
    d = length(func)
    upper_to_lower, lower_to_upper = _upper_to_lower_triangular_permutation(d)
    scale_factor = fill(scale, d)
    for i in 1:d
        if MOI.Utilities.is_diagonal_vectorized_index(i)
            scale_factor[i] = 1.0
        end
    end
    if moi_to_scs
        return (func.*scale_factor)[lower_to_upper]
    else
        return func[upper_to_lower] .* scale_factor
    end
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
