# Utilities for applying a scaled permutation to a `Vector`,
# `MOI.VectorOfVariables` or `MOI.VectorAffineFunction`.
# We might consider having a similar functionality in `MOI.Utilities`
# that is faster than `MOI.Utilities.scalarize` -> permute -> scale ->
# `MOI.Utilities.vectorize`.

function _scaled_permutation(
    func::Vector{T},
    scale_factor::Vector{T},
    permutation,
    inverse_permutation,
    scale_before::Bool,
) where {T}
    return if scale_before
        (func .* scale_factor)[permutation]
    else
        func[permutation] .* scale_factor
    end
end

function _scaled_permutation(
    func::MOI.VectorOfVariables,
    scale_factor::Vector{T},
    permutation,
    inverse_permutation,
    scale_before::Bool,
) where {T}
    n = length(func.variables)
    return MOI.VectorAffineFunction{T}(
        [MOI.VectorAffineTerm(
            i,
            MOI.ScalarAffineTerm(
                (scale_before ? scale_factor[permutation[i]] : scale_factor[i]),
                func.variables[permutation[i]],
            ),
        ) for i in 1:n],
        zeros(T, n),
    )
end

function _scaled_permutation(
    func::MOI.VectorAffineFunction{T},
    scale_factor::Vector{T},
    permutation,
    inverse_permutation,
    scale_before::Bool,
) where {T}
    scaled_constants = _scaled_permutation(
        func.constants,
        scale_factor,
        permutation,
        inverse_permutation,
        scale_before,
    )
    scaled_terms = map(func.terms) do term
        row = term.output_index
        i = inverse_permutation[row]
        scale_i = scale_before ? scale_factor[row] : scale_factor[i]
        return MOI.VectorAffineTerm(
            i,
            MOI.ScalarAffineTerm(
                term.scalar_term.coefficient * scale_i,
                term.scalar_term.variable,
            ),
        )
    end
    return MOI.VectorAffineFunction(scaled_terms, scaled_constants)
end

# TODO define `MOI.output_dimension` for `Vector`
_output_dimension(func::Vector) = length(func)
_output_dimension(func) = MOI.output_dimension(func)
