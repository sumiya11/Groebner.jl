
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
    graph=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(:linear_algebra!, matrix)

    rng = params.rng
    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end

    flag = if !isnothing(graph)
        linear_algebra!(graph, matrix, basis, linalg, arithmetic, rng)
    else
        linear_algebra!(matrix, basis, linalg, arithmetic, rng)
    end

    flag
end

function linear_algebra_autoreduce_basis!(
    matrix::MacaulayMatrix,
    basis::Basis,
    params::AlgorithmParameters,
    graph=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(:linear_algebra_autoreduce_basis!, matrix)

    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end

    flag = if !isnothing(graph)
        linear_algebra_autoreduce_basis!(graph, matrix, basis, linalg, arithmetic)
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
    @log level = -3 repr_matrix(matrix)

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
    @log level = -3 repr_matrix(matrix)

    reduce_matrix_lower_part_any_nonzero!(matrix, basis, arithmetic)
end

###

function linear_algebra!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :deterministic
        deterministic_sparse_linear_algebra!(matrix, basis, linalg, arithmetic)
    elseif linalg.algorithm === :randomized
        randomized_sparse_linear_algebra!(matrix, basis, linalg, arithmetic, rng)
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

# Linear algebra with a learned graph of computation
function linear_algebra!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :learn
        learn_sparse_linear_algebra!(graph, matrix, basis, arithmetic)
    else
        @assert linalg.algorithm === :apply
        apply_sparse_linear_algebra!(graph, matrix, basis, arithmetic)
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
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    reduce_matrix_lower_part!(matrix, basis, arithmetic)
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
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    randomized_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
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
    @log level = -3 repr_matrix(matrix)
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
    @log level = -3 repr_matrix(matrix)

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
    @log level = -3 repr_matrix(matrix)

    # Produce the RREF of AB
    interreduce_matrix_upper_part_sparsedense!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    reduce_matrix_lower_part_sparsedense!(matrix, basis, arithmetic)
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function learn_sparse_linear_algebra!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "learn_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    learn_reduce_matrix_lower_part!(graph, matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function apply_sparse_linear_algebra!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    # NOTE: here, we do not need to sort the rows in the upper part, as they
    # have already been collected in the right order
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "apply_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    flag = apply_reduce_matrix_lower_part!(graph, matrix, basis, arithmetic)
    if !flag
        return flag
    end
    # Interreduce CD
    apply_interreduce_matrix_pivots!(graph, matrix, basis, arithmetic)

    true
end

function linear_algebra_autoreduce_basis!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)
    deterministic_sparse_interreduction!(matrix, basis, arithmetic)
    true
end

function linear_algebra_autoreduce_basis!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)

    if linalg.algorithm === :learn
        learn_deterministic_sparse_interreduction!(graph, matrix, basis, arithmetic)
    else
        @assert linalg.algorithm === :apply
        apply_deterministic_sparse_interreduction!(graph, matrix, basis, arithmetic)
    end

    true
end

function deterministic_sparse_interreduction!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=true)

    true
end

function learn_deterministic_sparse_interreduction!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "learn_deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    learn_interreduce_matrix_pivots!(graph, matrix, basis, arithmetic, reversed_rows=true)

    true
end

function apply_deterministic_sparse_interreduction!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log level = -3 "apply_deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    flag = apply_interreduce_matrix_pivots!(
        graph,
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

# Normalize the row to have the leading coefficient equal to 1
function normalize_row!(
    row::Vector{T},
    arithmetic::A,
    first_nnz_index::Int=1
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @invariant @inbounds !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    @inbounds if isone(lead)
        return row
    end

    @inbounds pinv = mod_x(invmod(lead, divisor(arithmetic)), arithmetic)
    @inbounds row[1] = one(T)
    @inbounds for i in 2:length(row)
        row[i] = mod_x(row[i] * pinv, arithmetic)
    end

    row
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_row!(
    row::Vector{T},
    arithmetic::A,
    first_nnz_index::Int=1
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @invariant @inbounds !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    @inbounds if isone(lead)
        return row
    end

    @inbounds pinv = inv(lead)
    @inbounds row[1] = one(T)
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end

    row
end

# Linear combination of dense vector and sparse vector
function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::A
) where {I, T <: CoeffFF, A <: AbstractArithmeticZp}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    @inbounds mul = divisor(arithmetic) - row[indices[1]]
    # On our benchmarks, usually,
    #   length(row) / length(indices)
    # roughly varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = mod_x(row[idx] + mul * coeffs[j], arithmetic)
    end

    sum(rand(100))
end

# Linear combination of dense vector and sparse vector
function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::A
) where {I, T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
    end

    nothing
end

# Linear combination of two dense vectors
function reduce_dense_row_by_dense_row!(
    row1::Vector{T},
    row2::Vector{T},
    mul::T,
    arithmetic::A
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @invariant length(row1) == length(row2)

    @inbounds mul = divisor(arithmetic) - mul

    @inbounds for j in 1:length(row1)
        row1[j] = mod_x(row1[j] + mul * row2[j], arithmetic)
    end

    nothing
end

# Linear combination of two dense vectors
function reduce_dense_row_by_dense_row!(
    row1::Vector{T},
    row2::Vector{T},
    mul::T,
    arithmetic::A
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @invariant length(row1) == length(row2)

    @inbounds for j in 1:length(row1)
        row1[j] = row1[j] - mul * row2[j]
    end

    nothing
end

# Load the coefficients from `coeffs` into dense `row` at `indices`. Zero the
# entries of `row` before that.
function load_sparse_row!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T}
) where {I, T <: CoeffFF}
    @invariant length(indices) == length(coeffs)

    @inbounds for i in 1:length(row)
        row[i] = T(0)
    end

    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end

    nothing
end

# Load the coefficients from `coeffs` into dense `row` at `indices`. Zero the
# entries of `row` before that.
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
    row::Vector{T},
    from::J,
    to::J
) where {I, J, T <: Coeff}
    # NOTE: assumes the provided sparse row has the necessary capacity allocated
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

# Given a matrix of form 
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, the matrix of the following form is produced:
#   A B
#   0 D' 
@timeit function reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    _, nlow = nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    # TODO: can be allocated in the matrix once and for all iterations
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
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
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    true
end

@timeit function reduce_matrix_lower_part_sparsedense!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    nleft, nright = ncols_filled(matrix)
    _, nlower = nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.D_coeffs_dense, nlow)

    # Allocate the buffers
    row1 = zeros(C, nleft)

    @inbounds for i in 1:nlower
        row2 = zeros(C, nright)
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
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nup, _ = nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.upper_coeffs, nup)
    resize!(matrix.some_coeffs, matrix.nrows_filled_lower)

    # Allocate the buffers
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)

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

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    matrix.upper_part_is_rref = true

    true
end

@timeit function interreduce_matrix_upper_part_sparsedense!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    nrows, ncols = size(matrix)
    nleft, nright = ncols_filled(matrix)
    nup, _ = nrows_filled(matrix)

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
    row1 = zeros(C, nleft)

    @inbounds for i in nup:-1:1
        row2 = zeros(C, nright)

        nnz_column_indices = matrix.upper_rows[i]
        nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]

        @invariant length(nnz_column_indices) == length(nnz_coeffs)
        @invariant isone(nnz_coeffs[1])

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
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nleft, nright = ncols_filled(matrix)
    nupper, _ = nrows_filled(matrix)

    # Prepare the matrix
    resize!(matrix.lower_rows, nright)
    pivots = matrix.pivots
    new_pivots = 0
    any_zeroed = false

    # Allocate the buffers
    row = zeros(C, ncols)
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

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
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

    true, any_zeroed, not_reduced_to_zero
end

@timeit function interreduce_matrix_pivots_dense!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nleft, nright = ncols_filled(matrix)
    _, nlower = nrows_filled(matrix)

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
            reduce_dense_row_by_dense_row!(to_be_reduced, row2, mul, arithmetic)
        end

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
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
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    _, nlow = nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)

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

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    matrix.npivots = matrix.nrows_filled_lower = matrix.nrows_filled_lower
end

function reduce_matrix_lower_part_any_nonzero!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    _, nlow = nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)

    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = matrix.some_coeffs[row_idx_to_coeffs[i]]

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
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A,
    rng::AbstractRNG
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    _, nlow = nrows_filled(matrix)

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
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    rows_multipliers = zeros(C, rowsperblock)

    for i in 1:nblocks
        nrowsupper = nlow > i * rowsperblock ? i * rowsperblock : nlow
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        new_pivots_count = 0
        @inbounds while new_pivots_count < nrowstotal
            # Produce a random linear combination of several rows from the lower
            # part of the matrix
            for j in 1:nrowstotal
                # TODO: does not work for the rationals
                rows_multipliers[j] = mod_x(rand(rng, C), arithmetic)
            end
            row .= C(0)
            first_nnz_col = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                nnz_column_indices = matrix.lower_rows[rowidx]
                nnz_coeffs = basis.coeffs[row_idx_to_coeffs[rowidx]]
                @invariant length(nnz_column_indices) == length(nnz_coeffs)

                first_nnz_col = min(first_nnz_col, nnz_column_indices[1])
                for l in 1:length(nnz_coeffs)
                    colidx = nnz_column_indices[l]
                    row[colidx] =
                        mod_x(row[colidx] + rows_multipliers[k] * nnz_coeffs[l], arithmetic)
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

            new_column_indices, new_coeffs = new_empty_sparse_row(C)
        end
    end

    true
end

function randomized_hashcolumns_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A,
    rng::AbstractRNG
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nup, nlow = nrows_filled(matrix)
    nleft, nright = ncols_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.B_coeffs_dense, ncols)

    # Compute hashes of rows in the B block
    hash_vector = matrix.buffer_hash_vector
    resize!(hash_vector, nright)
    @inbounds for i in 1:nright
        hash_vector[i] = mod_x(rand(C), arithmetic)
    end
    @inbounds for i in 1:nup
        row_hash = zero(C)
        nnz_column_indices = matrix.upper_rows[i]
        nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        # @info "" nnz_column_indices nnz_coeffs
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                continue
            end
            row_hash =
                mod_x(row_hash + nnz_coeffs[j] * hash_vector[idx - nleft], arithmetic)
        end
        matrix.B_coeffs_dense[nnz_column_indices[1]] = [row_hash]
    end

    # Allocate the buffers
    row = zeros(C, ncols)
    row1 = zeros(C, nleft)
    row2 = zeros(C, 1)

    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
        nnz_column_indices = matrix.lower_rows[i]
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        for j in 1:nleft
            row1[j] = zero(C)
        end
        row_hash = zero(C)
        for j in 1:length(nnz_column_indices)
            idx = nnz_column_indices[j]
            if idx <= nleft
                row1[idx] = nnz_coeffs[j]
            else
                row_hash =
                    mod_x(row_hash + nnz_coeffs[j] * hash_vector[idx - nleft], arithmetic)
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

            reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)

            reduce_dense_row_by_dense_row!(row2, reducer, mult, arithmetic)
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

        row_hash = zero(C)
        for j in 1:length(new_column_indices)
            idx = new_column_indices[j]
            if idx <= nleft
                continue
            end
            row_hash =
                mod_x(row_hash + new_coeffs[j] * hash_vector[idx - nleft], arithmetic)
        end
        matrix.B_coeffs_dense[new_column_indices[1]] = [row_hash]

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    true
end

###
# Learn & apply linear algebra

# Returns `false` if any row reduced to zero (since we expect that on the apply
# stage the rows are linearly independent)
function apply_reduce_matrix_lower_part!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    _, nlow = nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)

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

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    true
end

function learn_interreduce_matrix_pivots!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    # Perform interreduction
    flag, _, not_reduced_to_zero =
        interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=reversed_rows)
    !flag && return flag

    # Update the computation graph
    _, ncols = size(matrix)
    nup, nlow = nrows_filled(matrix)
    push!(
        graph.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (matrix.upper_to_coeffs[1:nup], matrix.upper_to_mult[1:nup])
    )
    push!(graph.matrix_lower_rows, (Vector{Int}(), Vector{Int}()))

    true
end

function apply_interreduce_matrix_pivots!(
    graph::ComputationGraphF4,
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
    nup, nlow = nrows_filled(matrix)

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

function absolute_index_pivots_in_interreduction!(matrix::MacaulayMatrix, basis::Basis)
    _, ncols = size(matrix)
    nup, nlow = nrows_filled(matrix)

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

# TODO: this function desperately needs a docstring!
function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{I},
    new_coeffs::Vector{C},
    row::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::A,
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
) where {I, C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nleft, _ = ncols_filled(matrix)

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

        reduce_dense_row_by_dense_row!(row2, reducer, mult, arithmetic)
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
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, ncols = size(matrix)
    nup, nlow = nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(C, ncols)
    not_reduced_to_zero = Vector{Int}()
    useful_reducers = Set{Tuple{Int, Int, MonomIdx}}()
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    reducer_rows = Tuple{Int, Int, MonomIdx}[]

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

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    # Update the tracer information
    # NOTE: we sort reducers by their original position in the array of pivots.
    # This way, the rows are already sorted at the apply stage.
    useful_reducers_sorted = sort(collect(useful_reducers), by=reducer -> reducer[1])
    push!(
        graph.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (map(f -> f[2], useful_reducers_sorted), map(f -> f[3], useful_reducers_sorted))
    )
    push!(
        graph.matrix_lower_rows,
        (row_idx_to_coeffs[not_reduced_to_zero], matrix.lower_to_mult[not_reduced_to_zero])
    )

    true
end

# function apply_sparse_rref!(
#     graph::ComputationGraphF4,
#     ring::PolyRing,
#     matrix::MacaulayMatrix{C},
#     basis::Basis{C},
#     rng::AbstractRNG
# ) where {C <: Coeff}
#     _, ncols = size(matrix)
#     nlow = matrix.nrows_filled_lower

#     arithmetic = select_arithmetic(matrix.some_coeffs, ring.ch)

#     # move known matrix pivots,
#     # no copy
#     pivs, l2c_tmp = get_absolute_pivots!(matrix)
#     @log level = -7 "absolute_pivots!" pivs l2c_tmp

#     rowidx2coef = matrix.lower_to_coeffs
#     matrix.lower_to_coeffs = l2c_tmp

#     # unknown pivots,
#     # we will modify them inplace when reducing by pivs
#     upivs = matrix.lower_rows
#     densecoeffs = zeros(C, ncols)

#     @inbounds for i in 1:nlow
#         # select next row to be reduced
#         # npiv ~ exponents
#         rowexps = upivs[i]

#         # corresponding coefficients from basis
#         # (no need to copy here)
#         cfsref = basis.coeffs[rowidx2coef[i]]

#         # we load coefficients into dense array
#         # into rowexps indices
#         load_sparse_row!(densecoeffs, rowexps, cfsref)

#         # reduce it with known pivots from matrix.upper_rows
#         # first nonzero in densecoeffs is at startcol position
#         startcol = rowexps[1]
#         zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
#             densecoeffs,
#             matrix,
#             basis,
#             pivs,
#             ColumnLabel(startcol),
#             ColumnLabel(ncols),
#             ColumnLabel(-1),
#             arithmetic
#         )
#         # if fully reduced
#         if zeroed
#             # then something has gone wrong in tracing 
#             return false
#         end

#         # matrix coeffs sparsely stores coefficients of new row
#         matrix.some_coeffs[i] = newcfs
#         # add new pivot at column index newrow[1]
#         #  (which is the first nnz column of newrow)
#         pivs[newrow[1]] = newrow
#         # set ref to coefficient to matrix
#         # guaranteed to be from lower part
#         matrix.lower_to_coeffs[newrow[1]] = i

#         # normalize if needed
#         normalize_row!(matrix.some_coeffs[i], arithmetic)
#     end

#     interreduce_lower_part!(matrix, basis, pivs, arithmetic)
#     true
# end

###
# Re-enumerating columns in the matrix and other auxiliaries
# TODO: probably move this out of here?

function column_to_monom_mapping!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)
    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    column_to_monom = Vector{MonomIdx}(undef, load - 1)
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

    resize_basis_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.nprocessed

    @inbounds for i in 1:(matrix.npivots)
        colidx = rows[i][1]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.column_to_monom)
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
    resize_basis_if_needed!(basis, matrix.npivots)

    @inbounds for i in 1:(matrix.npivots)
        basis.nprocessed += 1
        basis.nnonredundant += 1
        basis.nonredundant[basis.nnonredundant] = basis.nprocessed
        if isassigned(matrix.some_coeffs, i)
            row = matrix.lower_rows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[basis.nprocessed] = matrix.some_coeffs[i]
            basis.monoms[basis.nprocessed] = row
        else
            empty!(basis.coeffs[basis.nprocessed])
            empty!(basis.monoms[basis.nprocessed])
        end
    end

    nothing
end

function insert_in_basis_hash_table_pivots(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomIdx}
) where {M}
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

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask, sdata[hidx].deg)

        ht.load += 1
    end

    nothing
end
