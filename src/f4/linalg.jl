
"""
    linear_algebra!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

linearly reduces the `CD` part with respect to the `AB` part, and interreduces
the result inplace. 

In other words, computes the reduced row echelon form of

    D - C inv(A) B

Returns `true` if successful and `false` otherwise.
"""
function linear_algebra! end

"""
    linear_algebra_autoreduce_basis!

Interreduces the rows of the given `MacaulayMatrix`.

Returns `true` if successful and `false` otherwise.
"""
function linear_algebra_autoreduce_basis! end

"""
    linear_algebra_normalform!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

linearly reduces the `CD` part with respect to the `AB` part.

In contrast to `linear_algebra!`, does not perform the final interreduction of
the matrix rows.

Returns `true` if successful and `false` otherwise.
"""
function linear_algebra_normalform! end

"""
    linear_algebra_isgroebner!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

checks that

    D - C inv(A) B

is zero.
"""
function linear_algebra_isgroebner! end

@noinline __throw_linalg_error(msg) = throw(DomainError("Linear algebra error: $msg"))

###
# Linear algebra, high level, main entry points

function linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    params::AlgorithmParameters,
    trace=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(:linear_algebra!, matrix)

    @stat matrix_block_sizes = matrix_block_sizes(matrix)
    @stat matrix_nnz = sum(i -> length(matrix.upper_rows[i]), 1:length(matrix.upper_rows))

    rng = params.rng
    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end

    threaded = if params.threaded_f4 === :yes && nthreads() > 1
        :yes
    elseif params.threaded_multimodular === :yes &&
           linalg.algorithm === :learn &&
           nthreads() > 1
        :yes
    elseif params.threaded_f4 === :auto
        # Multi-threading is enabled by default!
        :yes
    else
        :no
    end

    @log level = -4 "Calling linear algebra with" linalg threaded arithmetic

    flag = if !isnothing(trace)
        linear_algebra_with_trace!(trace, matrix, basis, linalg, threaded, arithmetic, rng)
    else
        linear_algebra!(matrix, basis, linalg, threaded, arithmetic, rng)
    end

    flag
end

function linear_algebra_autoreduce_basis!(
    matrix::MacaulayMatrix,
    basis::Basis,
    params::AlgorithmParameters,
    trace=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(:linear_algebra_autoreduce_basis!, matrix)

    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end

    flag = if !isnothing(trace)
        linear_algebra_autoreduce_basis!(trace, matrix, basis, linalg, arithmetic)
    else
        linear_algebra_autoreduce_basis!(matrix, basis, linalg, arithmetic)
    end

    flag
end

function linear_algebra_normalform!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @invariant matrix_well_formed(:linear_algebra_normalform!, matrix)

    sort_matrix_upper_rows!(matrix)
    @log level = -3 "linear_algebra_normalform!"
    @log level = -3 matrix_string_repr(matrix)

    reduce_matrix_lower_part_invariant_pivots!(matrix, basis, arithmetic)
end

function linear_algebra_isgroebner!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @invariant matrix_well_formed(:linear_algebra_isgroebner!, matrix)

    sort_matrix_upper_rows!(matrix)
    sort_matrix_lower_rows!(matrix)
    @log level = -3 "linear_algebra_isgroebner!"
    @log level = -3 matrix_string_repr(matrix)

    reduce_matrix_lower_part_any_nonzero!(matrix, basis, arithmetic)
end

###
# Dispatching between linear algebra algorithms

function should_use_threading(matrix, linalg, threaded)
    @assert threaded !== :auto
    nup, nlow, nleft, nright = matrix_block_sizes(matrix)
    if threaded === :yes
        # Opt for single thread for small matrices
        if nlow < 2 * nthreads()
            false
        elseif nleft + nright < 1_000
            false
        else
            true
        end
    else
        false
    end
end

function linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    threaded::Symbol,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :deterministic
        if should_use_threading(matrix, linalg, threaded)
            deterministic_sparse_linear_algebra_threaded!(matrix, basis, linalg, arithmetic)
        else
            deterministic_sparse_linear_algebra!(matrix, basis, linalg, arithmetic)
        end
    elseif linalg.algorithm === :randomized
        if should_use_threading(matrix, linalg, threaded)
            randomized_sparse_linear_algebra_threaded!(matrix, basis, linalg, arithmetic, rng)
        else
            randomized_sparse_linear_algebra!(matrix, basis, linalg, arithmetic, rng)
        end
    elseif linalg.algorithm === :experimental_1
        if linalg.sparsity === :sparse
            direct_rref_sparse_linear_algebra!(matrix, basis, linalg, arithmetic)
        else
            direct_rref_sparsedense_linear_algebra!(matrix, basis, linalg, arithmetic)
        end
    elseif linalg.algorithm === :experimental_2
        randomized_hashcolumns_sparse_linear_algebra!(matrix, basis, linalg, arithmetic, rng)
    else
        __throw_linalg_error("Not recognized linear algebra option: $(linalg.algorithm)")
    end

    flag
end

# Linear algebra with a learned trace of computation
function linear_algebra_with_trace!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    threaded::Symbol,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :learn
        if false && threaded === :yes
            learn_sparse_linear_algebra_threaded!(trace, matrix, basis, arithmetic)
        else
            learn_sparse_linear_algebra!(trace, matrix, basis, arithmetic)
        end
    else
        @assert linalg.algorithm === :apply
        apply_sparse_linear_algebra!(trace, matrix, basis, arithmetic)
    end

    flag
end

function deterministic_sparse_linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "deterministic_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    reduce_matrix_lower_part!(matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function deterministic_sparse_linear_algebra_threaded!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "deterministic_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    if false
        reduce_matrix_lower_part_threaded_lock_free!(matrix, basis, arithmetic)
    else
        reduce_matrix_lower_part_threaded_cas!(matrix, basis, arithmetic)
    end
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function randomized_sparse_linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "randomized_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    randomized_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function randomized_sparse_linear_algebra_threaded!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "randomized_sparse_linear_algebra_threaded!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    if true
        randomized_reduce_matrix_lower_part_threaded_cas!(matrix, basis, arithmetic, rng)
    else
        randomized_reduce_matrix_lower_part_threaded_cas_fair!(
            matrix,
            basis,
            arithmetic,
            rng
        )
    end
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function randomized_hashcolumns_sparse_linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "randomized_hashcolumns_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    randomized_hashcolumns_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function direct_rref_sparse_linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log level = -3 "direct_rref_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)

    # Produce the RREF of AB
    interreduce_matrix_upper_part!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    reduce_matrix_lower_part!(matrix, basis, arithmetic)
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function direct_rref_sparsedense_linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log level = -3 "direct_rref_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)

    # Produce the RREF of AB
    interreduce_matrix_upper_part_sparsedense!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    reduce_matrix_lower_part_sparsedense!(matrix, basis, arithmetic)
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function learn_sparse_linear_algebra!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "learn_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    learn_reduce_matrix_lower_part!(trace, matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function learn_sparse_linear_algebra_threaded!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "learn_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    learn_reduce_matrix_lower_part_threaded!(trace, matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function apply_sparse_linear_algebra!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    # NOTE: here, we do not need to sort the rows in the upper part, as they
    # have already been collected in the right order
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "apply_sparse_linear_algebra!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    flag = apply_reduce_matrix_lower_part!(trace, matrix, basis, arithmetic)
    if !flag
        return flag
    end
    # Interreduce CD
    apply_interreduce_matrix_pivots!(trace, matrix, basis, arithmetic)

    true
end

function linear_algebra_autoreduce_basis!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)
    deterministic_sparse_interf4_reduction!(matrix, basis, arithmetic)
    true
end

function linear_algebra_autoreduce_basis!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)

    if linalg.algorithm === :learn
        learn_deterministic_sparse_interf4_reduction!(trace, matrix, basis, arithmetic)
    else
        @assert linalg.algorithm === :apply
        apply_deterministic_sparse_interf4_reduction!(trace, matrix, basis, arithmetic)
    end

    true
end

function deterministic_sparse_interf4_reduction!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "deterministic_sparse_interf4_reduction!"
    @log level = -3 matrix_string_repr(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interf4_reduction!(matrix, basis)
    # Interreduce AB
    interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=true)

    true
end

function learn_deterministic_sparse_interf4_reduction!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "learn_deterministic_sparse_interf4_reduction!"
    @log level = -3 matrix_string_repr(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interf4_reduction!(matrix, basis)
    # Interreduce AB
    learn_interreduce_matrix_pivots!(trace, matrix, basis, arithmetic, reversed_rows=true)

    true
end

function apply_deterministic_sparse_interf4_reduction!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "apply_deterministic_sparse_interf4_reduction!"
    @log level = -3 matrix_string_repr(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interf4_reduction!(matrix, basis)
    # Interreduce AB
    flag = apply_interreduce_matrix_pivots!(
        trace,
        matrix,
        basis,
        arithmetic,
        reversed_rows=true
    )

    flag
end

###
# Basic routines for dense and sparse vectors

function new_empty_sparse_row(::Type{C}) where {C <: Coeff}
    Vector{ColumnLabel}(), Vector{C}()
end

function dense_row_mod_p!(
    row::Vector{T},
    arithmetic::AbstractArithmeticZp{T},
    from::Int=1,
    to::Int=length(row)
) where {T <: Union{CoeffZp, CompositeCoeffZp}}
    # NOTE: This loop is usually successfully vectorized by the compiler for all
    # kinds of arithmetic. Still, the resulting native code may be not optimal
    @fastmath @inbounds for i in from:to
        row[i] = mod_p(row[i], arithmetic)
    end
    nothing
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_row!(
    row::Vector{T},
    arithmetic::AbstractArithmeticZp{A, T},
    first_nnz_index::Int=1
) where {A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant @inbounds !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    isone(lead) && return row

    @inbounds pinv = inv_mod_p(A(lead), arithmetic) % T
    @inbounds row[1] = one(T)
    @inbounds for i in 2:length(row)
        row[i] = mod_p(A(row[i]) * A(pinv), arithmetic) % T
    end

    row
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_row!(
    row::Vector{T},
    arithmetic::AbstractArithmeticQQ{T},
    first_nnz_index::Int=1
) where {T <: CoeffQQ}
    @invariant @inbounds !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    isone(lead) && return row

    @inbounds pinv = inv(lead)
    @inbounds row[1] = one(T)
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end

    row
end

# Linear combination of dense vector and sparse vector modulo a prime.
# The most generic version.
function reduce_dense_row_by_sparse_row_mod_p!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticZp
) where {I, A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    @inbounds mul = divisor(arithmetic) - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        # if A === T, then the type cast is a no-op
        a = row[idx] + A(mul) * A(coeffs[j])
        row[idx] = mod_p(a, arithmetic)
    end

    nothing
end

# Linear combination of dense vector and sparse vector.
# The most generic version.
function reduce_dense_row_by_sparse_row!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticZp
) where {I, A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    @inbounds mul = divisor(arithmetic) - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        # if A === T, then the type cast is a no-op
        a = row[idx] + A(mul) * A(coeffs[j])
        row[idx] = a
    end

    nothing
end

# Linear combination of dense vector and sparse vector.
# Specialization for SignedArithmeticZp.
function reduce_dense_row_by_sparse_row!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::SignedArithmeticZp{A, T}
) where {I, A <: CoeffZp, T <: CoeffZp}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    # NOTE: mul is guaranteed to be < typemax(T)
    p2 = arithmetic.p2
    @inbounds mul = row[indices[1]] % T

    @fastmath @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] - A(mul) * A(coeffs[j])
        a = a + ((a >> (8 * sizeof(A) - 1)) & p2)
        row[idx] = a
    end

    nothing
end

# Linear combination of dense vector and sparse vector.
# Specialization for SignedCompositeArithmeticZp.
function reduce_dense_row_by_sparse_row!(
    row::Vector{CompositeInt{N, A}},
    indices::Vector{I},
    coeffs::Vector{CompositeInt{N, T}},
    arithmetic::SignedCompositeArithmeticZp{CompositeInt{N, A}, CompositeInt{N, T}}
) where {I, A <: CoeffZp, T <: CoeffZp, N}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    # NOTE: mul is guaranteed to be < typemax(T)
    p2 = arithmetic.p2s
    @inbounds mul = row[indices[1]].data .% T

    @fastmath @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx].data .- A.(mul) .* A.(coeffs[j].data)
        a = a .+ ((a .>> (8 * sizeof(A) - 1)) .& p2.data)
        row[idx] = CompositeInt(a)
    end

    nothing
end

# Linear combination of dense vector and sparse vector.
# Specialization for AbstractArithmeticQQ.
function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticQQ{T}
) where {I, T <: CoeffQQ}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
    end

    nothing
end

# Linear combination of two dense vectors modulo a prime
function reduce_dense_row_by_dense_row_mod_p!(
    row1::Vector{T},
    row2::Vector{T},
    mul::T,
    arithmetic::AbstractArithmeticZp
) where {T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant length(row1) == length(row2)

    @inbounds mul = divisor(arithmetic) - mul

    @inbounds for j in 1:length(row1)
        row1[j] = mod_p(row1[j] + mul * row2[j], arithmetic)
    end

    nothing
end

# Linear combination of two dense vectors
function reduce_dense_row_by_dense_row_mod_p!(
    row1::Vector{T},
    row2::Vector{T},
    mul::T,
    arithmetic::AbstractArithmeticQQ
) where {T <: CoeffQQ}
    @invariant length(row1) == length(row2)

    @inbounds for j in 1:length(row1)
        row1[j] = row1[j] - mul * row2[j]
    end

    nothing
end

# Load the contiguous coefficients from `coeffs` into dense `row` at `indices`.
# Zero the entries of `row` before that.
function load_sparse_row!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T}
) where {I, A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant length(indices) == length(coeffs)

    @inbounds for i in 1:length(row)
        row[i] = zero(A)
    end

    @inbounds for j in 1:length(indices)
        row[indices[j]] = A(coeffs[j])
    end

    nothing
end

# Load the contiguous coefficients from `coeffs` into dense `row` at `indices`.
# Zero the entries of `row` before that.
function load_sparse_row!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T}
) where {I, T <: CoeffQQ}
    @invariant length(indices) == length(coeffs)

    row .= T(0)
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end

    nothing
end

# Traverses the dense `row` at positions `from..to` and extracts all nonzero
# entries to the given sparse row. Returns the number of extracted nonzero
# elements.
function extract_sparse_row!(
    indices::Vector{I},
    coeffs::Vector{T},
    row::Vector{A},
    from::J,
    to::J
) where {I, J, T <: Coeff, A <: Coeff}
    # NOTE: also assumes that the provided sparse row has the necessary capacity
    # allocated
    @invariant length(indices) == length(coeffs)
    @invariant 1 <= from <= to <= length(row)

    z = zero(A)
    j = 1
    @inbounds for i in from:to
        if row[i] != z
            indices[j] = i
            # row[i] is expected to be less than typemax(T)
            coeffs[j] = row[i] % T
            j += 1
        end
    end

    j - 1
end

# Traverses the dense `row` at positions `from..to` and extracts all nonzero
# entries to the given sparse row. Returns the number of extracted nonzero
# elements.
# Specialization with eltype(coeffs) == eltype(row)
function extract_sparse_row!(
    indices::Vector{I},
    coeffs::Vector{T},
    row::Vector{T},
    from::J,
    to::J
) where {I, J, T <: Coeff}
    # NOTE: also assumes that the provided sparse row has the necessary capacity
    # allocated
    @invariant length(indices) == length(coeffs)
    @invariant 1 <= from <= to <= length(row)

    z = zero(T)
    j = 1
    @inbounds for i in from:to
        if row[i] != z
            indices[j] = i
            coeffs[j] = row[i]
            j += 1
        end
    end

    j - 1
end

# TODO: not really important atm
# function extract_sparse_row_2!(vals::Vector{UInt}, dense::Vector{UInt})
#     j = 1
#     z = Vec{4, UInt}((0, 0, 0, 0))
#     @inbounds for i in 1:8:length(dense)
#         v1 = SIMD.vload(Vec{4, UInt}, dense, i)
#         v2 = SIMD.vload(Vec{4, UInt}, dense, i + 4)
#         mask1 = !iszero(v1)
#         mask2 = !iszero(v2)
#         j1 = sum(mask1)
#         j2 = sum(mask2)
#         vstorec(v1, vals, j, mask1)
#         vstorec(v2, vals, j + j1, mask2)
#         j += j1 + j2
#     end
#     nothing
# end

###
# Linear algebra, low level

# Given a matrix of the following form, where A is in REF,
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, a matrix of the following form is produced:
#   A B
#   0 D' 
@timeit function reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    # TODO: can be allocated in the matrix once and for all iterations of F4
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        # Select a row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Load the coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Reduce the row with respect to the known pivots from the upper part of
        # the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add the pivot to the pool of available pivots
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        zeroed && continue

        @invariant length(new_column_indices) == length(new_coeffs)
        normalize_row!(new_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    true
end

function _reduce_matrix_lower_part_threaded_cas_worker!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    row_index_to_coeffs,
    block_start::Int,
    block_end::Int
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    pivots = matrix.pivots

    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    sentinels = matrix.sentinels

    @inbounds for i in block_start:block_end
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Load the coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        success = false
        while !success

            # Reduce the row with respect to the known `pivots` from the upper part
            # of the matrix.
            # NOTE: this also does partial interreduction of the lower matrix rows.
            # Upon discovery of a new pivot from the lower part of the matrix, we
            # add the pivot to the pool of available pivots
            first_nnz_column = nnz_column_indices[1]
            zeroed = reduce_dense_row_by_pivots_sparse!(
                new_column_indices,
                new_coeffs,
                row,
                matrix,
                basis,
                pivots,
                first_nnz_column,
                ncols,
                arithmetic,
                tmp_pos=-1
            )

            # If the row is fully reduced
            zeroed && break

            @invariant length(new_coeffs) == length(new_column_indices)

            # Set a reference to the coefficients of this row in the matrix
            old, success = Atomix.replace!(
                Atomix.IndexableRef(sentinels, (Int(new_column_indices[1]),)),
                Int8(0),
                Int8(1),
                Atomix.seq_cst,
                Atomix.monotonic
            )

            if success
                @invariant iszero(old)
                normalize_row!(new_coeffs, arithmetic)
                matrix.some_coeffs[i] = new_coeffs
                matrix.lower_to_coeffs[new_column_indices[1]] = i
                pivots[new_column_indices[1]] = new_column_indices
            else
                nnz_column_indices = new_column_indices
                nnz_coeffs = new_coeffs
            end
        end

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end
end

@timeit function reduce_matrix_lower_part_threaded_cas!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Calculate the size of the block
    nblocks = nthreads()
    nblocks = min(nblocks, nlow)
    # nblocks = 1
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    @log level = -3 "" nblocks rem rowsperblock nlow

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.sentinels, ncols)
    @inbounds for i in 1:ncols
        matrix.sentinels[i] = 0
    end
    for i in 1:nup
        matrix.sentinels[matrix.upper_rows[i][1]] = 1
    end

    tasks = Vector{Task}(undef, nblocks)

    for i in 2:nblocks
        block_start = 1 + (i - 1) * rowsperblock
        block_end = min(nlow, i * rowsperblock)
        block_start > nlow && continue
        @invariant 1 <= block_start <= block_end <= nlow
        tasks[i] = Base.Threads.@spawn _reduce_matrix_lower_part_threaded_cas_worker!(
            matrix,
            basis,
            arithmetic,
            row_index_to_coeffs,
            block_start,
            block_end
        )
    end

    _reduce_matrix_lower_part_threaded_cas_worker!(
        matrix,
        basis,
        arithmetic,
        row_index_to_coeffs,
        1,
        rowsperblock
    )

    for i in 1:nblocks
        if isassigned(tasks, i)
            fetch(tasks[i])
        end
    end

    true
end

function _reduce_matrix_lower_part_threaded_lock_free_worker!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{CoeffType, AccumType},
    row_index_to_coeffs,
    block_start::Int,
    block_end::Int
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    pivots = matrix.pivots

    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in block_start:block_end
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Load the coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add the pivot to the pool of available pivots
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        zeroed && continue

        @invariant length(new_coeffs) == length(new_column_indices)
        normalize_row!(new_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_coeffs
        matrix.lower_rows[i] = new_column_indices
        matrix.sentinels[i] = 1

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end
end

@timeit function reduce_matrix_lower_part_threaded_lock_free!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Calculate the size of the block in threading
    nblocks = nthreads()
    nblocks = min(nblocks, nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    @log level = -3 "" nblocks rem rowsperblock nlow

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.sentinels, nlow)
    @inbounds for i in 1:nlow
        matrix.sentinels[i] = 0
    end

    tasks = Vector{Task}(undef, nblocks)

    for i in 2:nblocks
        block_start = 1 + (i - 1) * rowsperblock
        block_end = min(nlow, i * rowsperblock)
        block_start > nlow && continue
        @invariant 1 <= block_start <= block_end <= nlow
        tasks[i] = Base.Threads.@spawn _reduce_matrix_lower_part_threaded_lock_free_worker!(
            matrix,
            basis,
            arithmetic,
            row_index_to_coeffs,
            block_start,
            block_end
        )
    end

    _reduce_matrix_lower_part_threaded_lock_free_worker!(
        matrix,
        basis,
        arithmetic,
        row_index_to_coeffs,
        1,
        rowsperblock
    )

    for i in 1:nblocks
        if isassigned(tasks, i)
            fetch(tasks[i])
        end
    end

    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        if iszero(matrix.sentinels[i])
            continue
        end

        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = matrix.some_coeffs[i]
        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Load the coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add the pivot to the pool of available pivots
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        zeroed && continue

        @invariant length(new_coeffs) == length(new_column_indices)
        normalize_row!(new_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    true
end

@timeit function reduce_matrix_lower_part_sparsedense!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    nleft, nright = matrix_ncols_filled(matrix)
    _, nlower = matrix_nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.D_coeffs_dense, nlow)

    # Allocate the buffers
    row1 = zeros(AccumType, nleft)

    @inbounds for i in 1:nlower
        row2 = zeros(CoeffType, nright)
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        nnz_coeffs = basis.coeffs[matrix.lower_to_coeffs[i]]
        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Extract the sparse row into two dense arrays, so that the row equals
        # [row1..., row2...]
        for j in 1:nleft
            row1[j] = zero(C)
        end
        for j in 1:nright
            row2[j] = zero(C)
        end
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                row1[idx] = nnz_coeffs[j]
            else
                row2[idx - nleft] = nnz_coeffs[j]
            end
        end

        # Reduce the rows simultaneously
        first_nnz_column = nnz_column_indices[1]
        pivot = reduce_dense_row_by_pivots_sparsedense!(
            row1,
            row2,
            matrix,
            basis,
            matrix.pivots,
            first_nnz_column,
            nleft,
            arithmetic
        )

        # If the row is fully reduced
        iszero(pivot) && continue

        normalize_row!(row2, arithmetic, pivot)
        matrix.D_coeffs_dense[i] = row2
        matrix.pivot_indices[i] = pivot
    end

    true
end

# Puts the AB part of the matrix into the RREF inplace.
@timeit function interreduce_matrix_upper_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, _ = matrix_nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.upper_coeffs, nup)
    resize!(matrix.some_coeffs, matrix.nrows_filled_lower)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in nup:-1:1
        # Locate the support and the coefficients of the row from the upper part
        # of the matrix
        nnz_column_indices = matrix.upper_rows[i]
        nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]

        # Extract the coefficients into a dense array
        @invariant isone(nnz_coeffs[1])
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: note the `tmp_pos=first_nnz_column` argument. It ensures that we
        # do not reduce the row with itself. Maybe think of a better name?..
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            matrix.upper_rows,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=first_nnz_column,
            computing_rref=true
        )

        # If the row is fully reduced
        if zeroed
            return false
        end

        @invariant length(new_coeffs) == length(new_column_indices)
        normalize_row!(new_coeffs, arithmetic)

        # Update the support and the coefficients of the vector
        matrix.upper_coeffs[i] = new_coeffs
        matrix.upper_rows[i] = new_column_indices

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    matrix.upper_part_is_rref = true

    true
end

@timeit function interreduce_matrix_upper_part_sparsedense!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    nrows, ncols = size(matrix)
    nleft, nright = matrix_ncols_filled(matrix)
    nup, _ = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots = Vector{Vector{ColumnLabel}}(undef, ncols)
    @inbounds for i in 1:nup
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
    end
    matrix.pivots = pivots

    resize!(matrix.pivot_indices, nrows)
    @inbounds for i in 1:nrows
        matrix.pivot_indices[i] = 0
    end

    resize!(matrix.upper_coeffs, nup)
    resize!(matrix.B_coeffs_dense, ncols)
    resize!(matrix.some_coeffs, matrix.nrows_filled_lower)

    # Allocate the buffers
    row1 = zeros(AccumType, nleft)

    @inbounds for i in nup:-1:1
        row2 = zeros(AccumType, nright)

        nnz_column_indices = matrix.upper_rows[i]
        nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]

        @invariant length(nnz_column_indices) == length(nnz_coeffs)
        @invariant isone(nnz_coeffs[1])

        # Extract the sparse row into two dense arrays, so that the row equals
        # [row1..., row2...]
        for j in 1:nleft
            row1[j] = zero(AccumType)
        end
        for j in 1:nright
            row2[j] = zero(AccumType)
        end
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                row1[idx] = nnz_coeffs[j]
            else
                row2[idx - nleft] = nnz_coeffs[j]
            end
        end

        first_nnz_column = nnz_column_indices[1]
        _ = reduce_dense_row_by_pivots_sparsedense!(
            row1,
            row2,
            matrix,
            basis,
            pivots,
            first_nnz_column + 1,
            nleft,
            arithmetic
        )

        matrix.B_coeffs_dense[first_nnz_column] = row2
    end

    matrix.upper_part_is_rref = true

    true
end

# Assumes that the pivots in the matrix are filled. Interreduces the pivots.
# Stores results into the lower part of the matrix.
#
# Returns the indices of rows (in the original coordinate system) that did not
# reduce to zero.
@timeit function interreduce_matrix_pivots!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType};
    reversed_rows::Bool=false
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nleft, nright = matrix_ncols_filled(matrix)
    nupper, _ = matrix_nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.lower_rows, nright)
    pivots = matrix.pivots
    new_pivots = 0
    any_zeroed = false

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    # Indices of rows that did no reduce to zero
    not_reduced_to_zero = Vector{Int}(undef, nright)

    # for each column in the block D..
    @inbounds for i in 1:nright
        abs_column_idx = ncols - i + 1
        # Check if there is a row that starts at `abs_column_idx`
        !isassigned(pivots, abs_column_idx) && continue

        # Locate the row support and coefficients
        nnz_column_indices = pivots[abs_column_idx]
        if abs_column_idx <= nleft
            # upper part of matrix
            nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[abs_column_idx]]
        else
            # lower part of matrix
            nnz_coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]]
        end
        @assert length(nnz_column_indices) == length(nnz_coeffs)

        # Load the row into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=first_nnz_column
        )

        # If the row reduced to zero
        if zeroed
            any_zeroed = true
            continue
        end

        new_pivots += 1
        not_reduced_to_zero[new_pivots] = i

        # Update the row support and coefficients.
        # TODO: maybe get rid of the reversed_rows parameter?
        if !reversed_rows
            matrix.lower_rows[new_pivots] = new_column_indices
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[new_pivots]
        else
            matrix.lower_rows[nupper - new_pivots + 1] = new_column_indices
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[nupper - new_pivots + 1]
        end
    end

    matrix.npivots = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    @log level = -3 "After interreduction" new_pivots

    true, any_zeroed, not_reduced_to_zero
end

@timeit function interreduce_matrix_pivots_dense!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType};
    reversed_rows::Bool=false
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nleft, nright = matrix_ncols_filled(matrix)
    _, nlower = matrix_nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.lower_to_coeffs, ncols)
    resize!(matrix.some_coeffs, ncols)
    resize!(matrix.lower_rows, ncols)
    pivot_indices = matrix.pivot_indices
    new_pivots = 0
    any_zeroed = false

    # Allocate the buffers
    not_reduced_to_zero = Vector{Int}(undef, nlower)

    @inbounds for i in 1:nlower
        if iszero(pivot_indices[i])
            any_zeroed = true
            continue
        end

        new_pivot = 0
        row2 = matrix.D_coeffs_dense[i]
        for j in 1:length(row2)
            if !iszero(row2[j])
                new_pivot = j
                break
            end
        end
        pivot_indices[i] = new_pivot

        if iszero(pivot_indices[i])
            any_zeroed = true
            continue
        end

        first_nnz_column = pivot_indices[i]
        normalize_row!(row2, arithmetic, first_nnz_column)

        for j in 1:nlower
            if j == i
                continue
            end
            if iszero(pivot_indices[j]) || pivot_indices[j] > first_nnz_column
                continue
            end
            to_be_reduced = matrix.D_coeffs_dense[j]
            if iszero(to_be_reduced[first_nnz_column])
                continue
            end
            mul = to_be_reduced[first_nnz_column]
            reduce_dense_row_by_dense_row_mod_p!(to_be_reduced, row2, mul, arithmetic)
        end

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
        resize!(new_column_indices, length(row2))
        resize!(new_coeffs, length(row2))
        cnt = extract_sparse_row!(
            new_column_indices,
            new_coeffs,
            row2,
            first_nnz_column,
            nright
        )
        resize!(new_column_indices, cnt)
        for j in 1:length(new_column_indices)
            new_column_indices[j] = nleft + new_column_indices[j]
        end
        resize!(new_coeffs, cnt)

        first_nnz_column = first_nnz_column + nleft
        @invariant isone(new_coeffs[1])

        new_pivots += 1
        not_reduced_to_zero[new_pivots] = first_nnz_column
        # update row and coeffs
        if !reversed_rows
            matrix.lower_rows[new_pivots] = new_column_indices
            matrix.lower_to_coeffs[first_nnz_column] = new_pivots
            matrix.some_coeffs[new_pivots] = new_coeffs
        else
            matrix.lower_rows[nlower - new_pivots + 1] = new_column_indices
            matrix.lower_to_coeffs[first_nnz_column] = new_pivots
            matrix.some_coeffs[new_pivots] = new_coeffs
        end
    end

    matrix.npivots = matrix.nrows_filled_lower = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    true, any_zeroed, not_reduced_to_zero
end

function reduce_matrix_lower_part_invariant_pivots!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = matrix.some_coeffs[row_idx_to_coeffs[i]]

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        first_nnz_column = nnz_column_indices[1]
        _ = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            # first_nnz_column,
            # TODO: using first_nnz_column is incorrect here; 
            # for the counter-example see
            # https://github.com/sumiya11/Groebner.jl/issues/82
            1,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # NOTE: no normalization here
        # NOTE: the array may be empty, signalling the row is zero
        @invariant length(new_coeffs) == length(new_column_indices)

        matrix.some_coeffs[i] = new_coeffs
        matrix.lower_rows[i] = new_column_indices
        matrix.lower_to_coeffs[i] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    matrix.npivots = matrix.nrows_filled_lower = matrix.nrows_filled_lower
end

function reduce_matrix_lower_part_any_nonzero!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If fully reduced
        zeroed && continue

        return false
    end

    true
end

###
# Randomized linear algebra

# The number of blocks to split the matrix into
function nblocks_in_randomized(nrows::Int)
    floor(Int, sqrt(nrows / 3)) + 1
end

# Given a matrix of form 
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, the matrix of the following form is produced:
#   A B
#   0 D' 
function randomized_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, matrix.nrows_filled_lower)

    # Set up the sizes of blocks
    nblocks = nblocks_in_randomized(nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem
    @log level = -3 """
    Rows in the lower part: $nlow
    The bumber of blocks: $nblocks
    Rows per block: $rowsperblock"""

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    rows_multipliers = zeros(AccumType, rowsperblock)

    for i in 1:nblocks
        nrowsupper = nlow > i * rowsperblock ? i * rowsperblock : nlow
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        new_pivots_count = 0
        @inbounds while new_pivots_count < nrowstotal
            # Produce a random linear combination of several rows from the lower
            # part of the matrix
            for j in 1:nrowstotal
                rows_multipliers[j] = mod_p(rand(rng, AccumType), arithmetic)
            end
            row .= AccumType(0)
            first_nnz_col = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                nnz_column_indices = matrix.lower_rows[rowidx]
                nnz_coeffs = basis.coeffs[row_idx_to_coeffs[rowidx]]
                @invariant length(nnz_column_indices) == length(nnz_coeffs)

                first_nnz_col = min(first_nnz_col, nnz_column_indices[1])
                for l in 1:length(nnz_coeffs)
                    colidx = nnz_column_indices[l]
                    # TODO: does not work for the rationals
                    # TODO: can vectorize random number generation
                    # TODO: consider sampling from a smaller range of values
                    row[colidx] = mod_p(
                        row[colidx] +
                        AccumType(rows_multipliers[k]) * AccumType(nnz_coeffs[l]),
                        arithmetic
                    )
                end
            end

            # Reduce the combination by rows from the upper part of the matrix
            zeroed = reduce_dense_row_by_pivots_sparse!(
                new_column_indices,
                new_coeffs,
                row,
                matrix,
                basis,
                pivots,
                first_nnz_col,
                ncols,
                arithmetic,
                tmp_pos=-1
            )

            zeroed && break
            new_pivots_count += 1

            normalize_row!(new_coeffs, arithmetic)
            @invariant length(new_column_indices) == length(new_coeffs)

            absolute_row_index = (i - 1) * rowsperblock + new_pivots_count
            matrix.some_coeffs[absolute_row_index] = new_coeffs
            pivots[new_column_indices[1]] = new_column_indices
            matrix.lower_to_coeffs[new_column_indices[1]] = absolute_row_index

            new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
        end
    end

    true
end

function randomized_reduce_matrix_lower_part_threaded_cas!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    nblocks = nblocks_in_randomized(nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.sentinels, ncols)
    sentinels = matrix.sentinels
    @inbounds for i in 1:ncols
        sentinels[i] = 0
    end
    for i in 1:nup
        sentinels[matrix.upper_rows[i][1]] = 1
    end

    rows_multipliers_buffers = map(_ -> zeros(AccumType, rowsperblock), 1:nthreads())
    row_buffers = map(_ -> zeros(AccumType, ncols), 1:nthreads())
    # thread_local_rngs = map(_ -> copy(rng), 1:nthreads())

    @inbounds Base.Threads.@threads :static for i in 1:nblocks
        nrowsupper = min(i * rowsperblock, nlow)
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        t_id = threadid()
        rows_multipliers = rows_multipliers_buffers[t_id]
        row = row_buffers[t_id]
        # new_column_indices, new_coeffs = new_row_buffers[t_id]
        # rng = thread_local_rngs[t_id]
        # row = zeros(AccumType, ncols)
        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
        # rows_multipliers = zeros(AccumType, rowsperblock)
        local_rng = copy(rng)

        new_pivots_count = 0
        @inbounds while new_pivots_count < nrowstotal
            # Produce a random linear combination of several rows from the lower
            # part of the matrix
            for j in 1:nrowstotal
                rows_multipliers[j] = mod_p(rand(local_rng, AccumType), arithmetic)
            end
            row .= AccumType(0)
            first_nnz_col = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                nnz_column_indices = matrix.lower_rows[rowidx]
                nnz_coeffs = basis.coeffs[row_index_to_coeffs[rowidx]]
                @invariant length(nnz_column_indices) == length(nnz_coeffs)

                first_nnz_col = min(first_nnz_col, nnz_column_indices[1])
                for l in 1:length(nnz_coeffs)
                    colidx = nnz_column_indices[l]
                    row[colidx] = mod_p(
                        row[colidx] +
                        AccumType(rows_multipliers[k]) * AccumType(nnz_coeffs[l]),
                        arithmetic
                    )
                end
            end

            success = false
            while !success
                # Reduce the combination by rows from the upper part of the matrix
                zeroed = reduce_dense_row_by_pivots_sparse!(
                    new_column_indices,
                    new_coeffs,
                    row,
                    matrix,
                    basis,
                    pivots,
                    first_nnz_col,
                    ncols,
                    arithmetic,
                    tmp_pos=-1
                )

                if zeroed
                    new_pivots_count = nrowstotal
                    break
                end

                @invariant length(new_column_indices) == length(new_coeffs)

                old, success = Atomix.replace!(
                    Atomix.IndexableRef(sentinels, (Int(new_column_indices[1]),)),
                    Int8(0),
                    Int8(1),
                    Atomix.acquire_release,
                    Atomix.monotonic
                )

                if success
                    @invariant iszero(old)

                    new_pivots_count += 1
                    absolute_row_index = (i - 1) * rowsperblock + new_pivots_count

                    normalize_row!(new_coeffs, arithmetic)
                    matrix.some_coeffs[absolute_row_index] = new_coeffs
                    matrix.lower_to_coeffs[new_column_indices[1]] = absolute_row_index
                    pivots[new_column_indices[1]] = new_column_indices
                else
                    first_nnz_col = new_column_indices[1]
                end
            end

            new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
        end
    end

    true
end

function randomized_hashcolumns_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)
    nleft, nright = matrix_ncols_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.B_coeffs_dense, ncols)

    # Compute hashes of rows in the B block
    hash_vector = matrix.buffer_hash_vector
    resize!(hash_vector, nright)
    @inbounds for i in 1:nright
        hash_vector[i] = mod_p(rand(AccumType), arithmetic)
    end
    @inbounds for i in 1:nup
        row_hash = zero(AccumType)
        nnz_column_indices = matrix.upper_rows[i]
        nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash + AccumType(nnz_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                arithmetic
            )
        end
        matrix.B_coeffs_dense[nnz_column_indices[1]] = [row_hash % CoeffType]
    end

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    row1 = zeros(AccumType, nleft)
    row2 = zeros(AccumType, 1)

    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        for j in 1:nleft
            row1[j] = zero(AccumType)
        end
        row_hash = zero(AccumType)
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                row1[idx] = nnz_coeffs[j]
            else
                row_hash = mod_p(
                    row_hash +
                    AccumType(nnz_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                    arithmetic
                )
            end
        end
        row2[1] = row_hash

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        first_nnz_column = nnz_column_indices[1]
        pivot = 0
        for q in first_nnz_column:ncols
            if iszero(row[q])
                continue
            end

            if !isassigned(pivots, q)
                if iszero(pivot)
                    pivot = q
                end
                continue
            end

            indices = pivots[q]
            if q <= nleft
                # if reducer is from the upper part of the matrix
                coeffs = basis.coeffs[matrix.upper_to_coeffs[q]]
            else
                # if reducer is from the lower part of the matrix
                coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[q]]
            end

            mult = row[q]
            reducer = matrix.B_coeffs_dense[q]

            reduce_dense_row_by_sparse_row_mod_p!(row, indices, coeffs, arithmetic)

            reduce_dense_row_by_dense_row_mod_p!(row2, reducer, mult, arithmetic)
        end
        zeroed = iszero(row2[1])

        # if fully reduced
        if zeroed && iszero(pivot)
            continue
        end

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # otherwise, reduce once again
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        if zeroed
            continue
        end

        normalize_row!(new_coeffs, arithmetic)

        matrix.some_coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        row_hash = zero(AccumType)
        for j in 1:length(new_column_indices)
            idx = new_column_indices[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash + AccumType(new_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                arithmetic
            )
        end
        matrix.B_coeffs_dense[new_column_indices[1]] = [row_hash]

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    true
end

###
# Learn & apply linear algebra

# Returns `false` if any row reduced to zero (since we expect that on the apply
# stage the rows are linearly independent)
function apply_reduce_matrix_lower_part!(
    trace::TraceF4,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]

        @invariant length(nnz_column_indices) == length(nnz_coeffs)

        # Load coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        if zeroed
            return false
        end

        @invariant length(new_coeffs) == length(new_column_indices)
        normalize_row!(new_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    true
end

function learn_interreduce_matrix_pivots!(
    trace::TraceF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    # Perform interreduction
    flag, _, not_reduced_to_zero =
        interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=reversed_rows)
    !flag && return flag

    # Update the computation trace
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)
    push!(
        trace.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    push!(trace.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        trace.matrix_upper_rows,
        (matrix.upper_to_coeffs[1:nup], matrix.upper_to_mult[1:nup])
    )
    push!(trace.matrix_lower_rows, (Vector{Int}(), Vector{Int}()))

    true
end

function apply_interreduce_matrix_pivots!(
    trace::TraceF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    flag, any_zeroed, _ =
        interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=reversed_rows)
    flag && !any_zeroed
end

###
# Utilities

function absolute_index_pivots!(matrix::MacaulayMatrix)
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    pivots = Vector{Vector{ColumnLabel}}(undef, ncols)
    @inbounds for i in 1:nup
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
    end

    absolute_lower_to_coeffs = Vector{ColumnLabel}(undef, max(ncols, nlow))
    @inbounds for i in 1:nlow
        absolute_lower_to_coeffs[matrix.lower_rows[i][1]] = matrix.lower_to_coeffs[i]
    end

    lower_to_coeffs = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = absolute_lower_to_coeffs
    matrix.pivots = pivots

    pivots, lower_to_coeffs
end

function absolute_index_pivots_in_interf4_reduction!(matrix::MacaulayMatrix, basis::Basis)
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    resize!(matrix.lower_rows, ncols)
    resize!(matrix.upper_to_coeffs, ncols)
    resize!(matrix.upper_to_mult, ncols)
    resize!(matrix.lower_to_coeffs, ncols)
    resize!(matrix.lower_to_mult, ncols)
    resize!(matrix.some_coeffs, ncols)

    pivots = Vector{Vector{ColumnLabel}}(undef, ncols)
    @inbounds for i in 1:(nup + nlow)
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
        matrix.lower_to_coeffs[matrix.upper_rows[i][1]] = i
        matrix.some_coeffs[i] = copy(basis.coeffs[matrix.upper_to_coeffs[i]])
    end

    matrix.pivots = pivots

    pivots
end

function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{I},
    new_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmeticZp{A, C},
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
) where {I, C <: Union{CoeffZp, CompositeCoeffZp}, A <: Union{CoeffZp, CompositeCoeffZp}}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero element in the reduced row
    new_pivot = -1

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.some_coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            if matrix.upper_part_is_rref || computing_rref
                coeffs = matrix.upper_coeffs[i]
            else
                coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            end
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        reduce_dense_row_by_sparse_row_mod_p!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_column_indices, nonzeros)
    resize!(new_coeffs, nonzeros)
    extract_sparse_row!(
        new_column_indices,
        new_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{I},
    new_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::DelayedArithmeticZp{A, C},
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
) where {I, C <: CoeffZp, A <: CoeffZp}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero element in the reduced row
    new_pivot = -1

    n_adds = 0
    n_safe = n_safe_consecutive_additions(arithmetic)

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        row[i] = mod_p(row[i], arithmetic)
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.some_coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            if matrix.upper_part_is_rref || computing_rref
                coeffs = matrix.upper_coeffs[i]
            else
                coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            end
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)

        n_adds += 1
        if iszero(n_adds % n_safe)
            dense_row_mod_p!(row, arithmetic, i, end_column)
        end
        row[i] = zero(row[i])

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_column_indices, nonzeros)
    resize!(new_coeffs, nonzeros)
    extract_sparse_row!(
        new_column_indices,
        new_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{I},
    new_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::Union{SignedArithmeticZp{A, C}, SignedCompositeArithmeticZp{A, C}},
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
) where {I, C <: Union{CoeffZp, CompositeCoeffZp}, A <: Union{CoeffZp, CompositeCoeffZp}}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero element in the reduced row
    new_pivot = -1

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end
        row[i] = mod_p(row[i], arithmetic)
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.some_coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            if matrix.upper_part_is_rref || computing_rref
                coeffs = matrix.upper_coeffs[i]
            else
                coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            end
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if nonzeros == 0
        return true
    end

    if end_column != length(row)
        dense_row_mod_p!(row, arithmetic, end_columns + 1, length(row))
    end

    # form the resulting row in sparse format
    resize!(new_column_indices, nonzeros)
    resize!(new_coeffs, nonzeros)
    extract_sparse_row!(
        new_column_indices,
        new_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{I},
    new_coeffs::Vector{C},
    row::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmeticQQ,
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
) where {I, C <: Coeff}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero element in the reduced row
    new_pivot = -1

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.some_coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            if matrix.upper_part_is_rref || computing_rref
                coeffs = matrix.upper_coeffs[i]
            else
                coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            end
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_column_indices, nonzeros)
    resize!(new_coeffs, nonzeros)
    extract_sparse_row!(
        new_column_indices,
        new_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

# TODO: this function desperately needs a docstring!
function reduce_dense_row_by_pivots_sparsedense!(
    row1::Vector{T},
    row2::Vector{T},
    matrix::MacaulayMatrix,
    basis::Basis,
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::A
) where {T <: Coeff, I, A <: AbstractArithmetic}
    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row1[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i)
            continue
        end

        mult = row1[i]
        reducer = matrix.B_coeffs_dense[i]

        reduce_dense_row_by_dense_row_mod_p!(row2, reducer, mult, arithmetic)
    end

    pivot = 0
    @inbounds for i in 1:length(row2)
        if !iszero(row2[i])
            pivot = i
            break
        end
    end

    pivot
end

###
# Linear algebra with tracer

record_active_reducer(active_reducers::Nothing, matrix, idx) = nothing
function record_active_reducer(active_reducers, matrix, idx)
    push!(active_reducers, (idx, matrix.upper_to_coeffs[idx], matrix.upper_to_mult[idx]))
    nothing
end

function learn_reduce_matrix_lower_part!(
    trace::TraceF4,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    not_reduced_to_zero = Vector{Int}()
    useful_reducers = Set{Tuple{Int, Int, MonomId}}()
    new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    reducer_rows = Tuple{Int, Int, MonomId}[]

    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        # Additionally record the indices of rows that participated in reduction
        # of the given row
        empty!(reducer_rows)
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            reducer_rows,
            tmp_pos=-1
        )

        # if fully reduced
        zeroed && continue

        # NOTE: we are not recording reducers from the lower part of the matrix
        push!(not_reduced_to_zero, i)
        for reducer_row in reducer_rows
            push!(useful_reducers, reducer_row)
        end

        @invariant length(new_column_indices) == length(new_coeffs)
        normalize_row!(new_coeffs, arithmetic)

        matrix.some_coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    # Update the tracer information
    # NOTE: we sort reducers by their original position in the array of pivots.
    # This way, the rows are already sorted at the apply stage.
    useful_reducers_sorted = sort(collect(useful_reducers), by=reducer -> reducer[1])
    push!(
        trace.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    push!(trace.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        trace.matrix_upper_rows,
        (map(f -> f[2], useful_reducers_sorted), map(f -> f[3], useful_reducers_sorted))
    )
    push!(
        trace.matrix_lower_rows,
        (row_idx_to_coeffs[not_reduced_to_zero], matrix.lower_to_mult[not_reduced_to_zero])
    )

    true
end

function learn_reduce_matrix_lower_part_threaded!(
    trace::TraceF4,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    if nthreads() == 1
        @log level = 0 """
        Using multi-threaded linear algebra with nthreads() == 1. 
        Something probably went wrong."""
    end

    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    resize!(matrix.sentinels, ncols)
    sentinels = matrix.sentinels
    @inbounds for i in 1:ncols
        sentinels[i] = 0
    end
    for i in 1:nup
        sentinels[matrix.upper_rows[i][1]] = 1
    end

    row_buffers = map(_ -> zeros(AccumType, ncols), 1:nthreads())
    useful_reducers_buffers = map(_ -> Set{Tuple{Int, Int, MonomId}}(), 1:nthreads())
    not_reduced_to_zero = zeros(Int8, nlow)

    @inbounds Base.Threads.@threads :static for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        t_id = threadid()
        row = row_buffers[t_id]
        useful_reducers = useful_reducers_buffers[t_id]
        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)

        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        success = false
        while !success
            first_nnz_col = nnz_column_indices[1]

            # Reduce the combination by rows from the upper part of the matrix
            zeroed = reduce_dense_row_by_pivots_sparse!(
                new_column_indices,
                new_coeffs,
                row,
                matrix,
                basis,
                pivots,
                first_nnz_col,
                ncols,
                arithmetic,
                useful_reducers,
                tmp_pos=-1
            )

            if zeroed
                break
            end

            @invariant length(new_column_indices) == length(new_coeffs)

            old, success = Atomix.replace!(
                Atomix.IndexableRef(sentinels, (Int(new_column_indices[1]),)),
                Int8(0),
                Int8(1),
                Atomix.acquire_release,
                Atomix.monotonic
            )

            if success
                @invariant iszero(old)

                normalize_row!(new_coeffs, arithmetic)
                matrix.some_coeffs[i] = new_coeffs
                matrix.lower_to_coeffs[new_column_indices[1]] = i
                pivots[new_column_indices[1]] = new_column_indices

                # NOTE: we are not recording reducers from the lower part of the matrix
                not_reduced_to_zero[i] = 1
            else
                nnz_column_indices = new_column_indices
                nnz_coeffs = new_coeffs
            end
        end

        new_column_indices, new_coeffs = new_empty_sparse_row(CoeffType)
    end

    # Update the tracer information
    # NOTE: we sort reducers by their original position in the array of pivots.
    # This way, the rows are already sorted at the apply stage.
    useful_reducers_sorted =
        sort(collect(reduce(union!, useful_reducers_buffers)), by=reducer -> reducer[1])
    push!(
        trace.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    not_reduced_to_zero_indices = findall(!iszero, not_reduced_to_zero)
    push!(trace.matrix_nonzeroed_rows, not_reduced_to_zero_indices)
    push!(
        trace.matrix_upper_rows,
        (map(f -> f[2], useful_reducers_sorted), map(f -> f[3], useful_reducers_sorted))
    )
    push!(
        trace.matrix_lower_rows,
        (
            row_idx_to_coeffs[not_reduced_to_zero_indices],
            matrix.lower_to_mult[not_reduced_to_zero_indices]
        )
    )

    true
end

###
# Re-enumerating columns in the matrix and other auxiliaries

@timeit function transform_polynomial_multiple_to_matrix_row!(
    matrix::MacaulayMatrix,
    symbolic_ht::MonomialHashtable{M},
    basis_ht::MonomialHashtable{M},
    htmp::MonomHash,
    etmp::M,
    poly::Vector{MonomId}
) where {M <: Monom}
    row = similar(poly)
    resize_hashtable_if_needed!(symbolic_ht, length(poly))

    insert_multiplied_poly_in_hashtable!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
end

function column_to_monom_mapping!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)
    @invariant !symbol_ht.frozen

    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    column_to_monom = Vector{MonomId}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        # column to hash index
        column_to_monom[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == PIVOT_COLUMN
            k += 1
        end
    end

    sort_columns_by_labels!(column_to_monom, symbol_ht)

    matrix.ncols_left = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.ncols_right = load - matrix.ncols_left - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        hv = hdata[column_to_monom[k]]
        hdata[column_to_monom[k]] = Hashvalue(k, hv.hash, hv.divmask, hv.deg)
    end

    @inbounds for k in 1:(matrix.nrows_filled_upper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nrows_filled_lower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.column_to_monom = column_to_monom
end

function convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    # We mutate the basis directly by adding new elements

    basis_resize_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.nprocessed

    @inbounds for i in 1:(matrix.npivots)
        colidx = rows[i][1]
        insert_in_basis_hashtable_pivots(rows[i], ht, symbol_ht, matrix.column_to_monom)
        basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
        basis.monoms[crs + i] = matrix.lower_rows[i]
        @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
    end

    basis.nfilled += matrix.npivots
end

function convert_rows_to_basis_elements_nf!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    basis_resize_if_needed!(basis, matrix.npivots)

    @inbounds for i in 1:(matrix.npivots)
        basis.nprocessed += 1
        basis.nnonredundant += 1
        basis.nonredundant[basis.nnonredundant] = basis.nprocessed
        if isassigned(matrix.some_coeffs, i)
            row = matrix.lower_rows[i]
            insert_in_basis_hashtable_pivots(row, ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[basis.nprocessed] = matrix.some_coeffs[i]
            basis.monoms[basis.nprocessed] = row
        else
            empty!(basis.coeffs[basis.nprocessed])
            empty!(basis.monoms[basis.nprocessed])
        end
    end

    nothing
end

function insert_in_basis_hashtable_pivots(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomId}
) where {M <: Monom}
    resize_hashtable_if_needed!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.monoms
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = column_to_monom[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = next_lookup_index(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if ishashcollision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        @invariant !ht.frozen

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask, sdata[hidx].deg)

        ht.load += 1
    end

    nothing
end