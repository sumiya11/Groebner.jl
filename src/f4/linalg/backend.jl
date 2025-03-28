# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_deterministic_sparse!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @invariant matrix.nrows_filled_upper == matrix.ncols_left
    for i in 1:(matrix.nrows_filled_upper)
        @invariant matrix.upper_rows[i][1] == i
    end

    # Reduce CD with AB
    linalg_reduce_matrix_lower_part!(matrix, basis, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function linalg_deterministic_sparse_interreduction!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)

    # Prepare the matrix
    linalg_prepare_matrix_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=true)
    true
end

###
# Low level

# Given a matrix of the following form, where A is in REF,
#   A B
#   C D
# reduces the lower part CD with respect to the pivots in the upper part AB.
# As a result, a matrix of the following form is produced:
#   A B
#   0 D'
# The new pivots in the D' block are known, but possibly not fully interreduced.
function linalg_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    # TODO: allocate just once, and then reuse for later iterations of F4
    row = [AccumType(zero(basis.coeffs[1][1])) for _ in 1:ncols]
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    # At each iteration, take a row from the lower part of the matrix (from the
    # CD block), reduce it by the pivots, and put the result back to the matrix
    @inbounds for i in 1:nlow
        # Select a row to be reduced from the lower part of the matrix
        # NOTE: no copy of coefficients is needed
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the coefficients into a dense array
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Reduce the row with respect to the known pivots
        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic
        )

        # If the row is fully reduced
        zeroed && continue

        @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)
        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        # Store the new row in the matrix, AND add the new row to the list of
        # known pivots
        pivots[new_sparse_row_support[1]] = new_sparse_row_support

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    end

    true
end

# Given a matrix of the form,
#   A B
#   0 D'
# interreduces the D' block. Stores results in the lower part of the matrix inplace.
#
# Assumes that the pivots in the matrix are correctly filled.
# 
# Returns the indices of the rows that did not reduce to zero.
function linalg_interreduce_matrix_pivots!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType};
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
    row = [AccumType(zero(basis.coeffs[1][1])) for _ in 1:ncols]
    # Indices of rows that did no reduce to zero
    not_reduced_to_zero = Vector{Int}(undef, nright)

    # for each column in the block D..
    @inbounds for i in 1:nright
        abs_column_idx = ncols - i + 1
        # Check if there is a row that starts at `abs_column_idx`
        !isassigned(pivots, abs_column_idx) && continue

        # Locate the row support and coefficients
        sparse_row_support = pivots[abs_column_idx]
        if abs_column_idx <= nleft
            # upper part of matrix
            sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[abs_column_idx]]
        else
            # lower part of matrix
            sparse_row_coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]]
        end
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the row into a dense array
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            ignore_column=first_nnz_column
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
            matrix.lower_rows[new_pivots] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[new_pivots]
        else
            matrix.lower_rows[nupper - new_pivots + 1] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[nupper - new_pivots + 1]
        end
    end

    matrix.npivots = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    true, any_zeroed, not_reduced_to_zero
end

function linalg_reduce_matrix_lower_part_do_not_modify_pivots!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = [AccumType(zero(basis.coeffs[1][1])) for _ in 1:ncols]
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = matrix.some_coeffs[row_idx_to_coeffs[i]]

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        first_nnz_column = sparse_row_support[1]
        _ = linalg_reduce_dense_row_by_pivots_sparse!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
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
            arithmetic
        )

        # NOTE: no normalization here
        # NOTE: the array may be empty, signalling the row is zero
        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        matrix.lower_rows[i] = new_sparse_row_support
        matrix.lower_to_coeffs[i] = i

        new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    end

    matrix.npivots = matrix.nrows_filled_lower = matrix.nrows_filled_lower
end

function linalg_reduce_matrix_lower_part_all_zero!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)

    # Allocate the buffers
    row = [AccumType(zero(basis.coeffs[1][1])) for _ in 1:ncols]
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic
        )

        # If fully reduced
        zeroed && continue

        return false
    end

    true
end

###
# Utilities

function linalg_prepare_matrix_pivots!(matrix::MacaulayMatrix)
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    pivots = Vector{Vector{ColumnLabel}}(undef, ncols)
    @inbounds for i in 1:nup
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
    end

    absolute_lower_to_coeffs = Vector{Int}(undef, max(ncols, nlow))
    @inbounds for i in 1:nlow
        absolute_lower_to_coeffs[matrix.lower_rows[i][1]] = matrix.lower_to_coeffs[i]
    end

    lower_to_coeffs = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = absolute_lower_to_coeffs
    matrix.pivots = pivots

    pivots, lower_to_coeffs
end

function linalg_prepare_matrix_pivots_in_interreduction!(matrix::MacaulayMatrix, basis::Basis)
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

###
# Reducing one dense vector by many pivots.

# Reduces the given dense row by the given pivots.
#
# Returns `true` if the row is reduced to zero and `false`, otherwise.
# Writes the nonzero elements of the resulting row to `new_sparse_row_support`
# and `new_sparse_row_coeffs`.
#
# The function may mutate
#   - `new_sparse_row_support`,
#   - `new_sparse_row_coeffs`, 
#   - `row`.
# It also mutates `active_reducers`, if provided.
function linalg_reduce_dense_row_by_pivots_sparse!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmeticZp{A, C},
    active_reducers=nothing;
    ignore_column::Integer=-1
) where {I, C <: Union{CoeffZp, CompositeCoeffZp}, A <: Union{CoeffZp, CompositeCoeffZp}}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    n_nonzeros = 0
    new_pivot_column = -1

    @inbounds for i in start_column:end_column
        # if the element is zero -- no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column equal to i
        if !isassigned(pivots, i) || ignore_column == i
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # At this point, we have determined that pivots[i] is a valid pivot.

        # Locate the support and the coefficients of the pivot row
        pivot_support = pivots[i]
        if i <= nleft
            # if pivot comes from the original pivots
            pivot_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            record_active_reducer(active_reducers, matrix, i)
        else
            pivot_coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(pivot_support) == length(pivot_coeffs)

        linalg_vector_addmul_sparsedense_mod_p!(row, pivot_support, pivot_coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all row elements reduced to zero!
    if n_nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_sparse_row_support, n_nonzeros)
    resize!(new_sparse_row_coeffs, n_nonzeros)
    linalg_extract_sparse_row!(
        new_sparse_row_support,
        new_sparse_row_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

# Specialization for delayed modular arithmetic
function linalg_reduce_dense_row_by_pivots_sparse!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::DelayedArithmeticZp{A, C},
    active_reducers=nothing;
    ignore_column::Integer=-1
) where {I, C <: CoeffZp, A <: CoeffZp}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    n_nonzeros = 0
    new_pivot_column = -1

    n_adds = 0
    n_safe = n_safe_consecutive_additions(arithmetic)

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        row[i] = mod_p(row[i], arithmetic)
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || ignore_column == i
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        linalg_vector_addmul_sparsedense!(row, indices, coeffs, arithmetic)

        n_adds += 1
        if iszero(n_adds % n_safe)
            linalg_dense_row_mod_p!(row, arithmetic, i, end_column)
        end
        row[i] = zero(row[i])

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if n_nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_sparse_row_support, n_nonzeros)
    resize!(new_sparse_row_coeffs, n_nonzeros)
    linalg_extract_sparse_row!(
        new_sparse_row_support,
        new_sparse_row_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

# Specialization for signed modular arithmetic
function linalg_reduce_dense_row_by_pivots_sparse!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::Union{SignedArithmeticZp{A, C}, SignedCompositeArithmeticZp{A, C}},
    active_reducers=nothing;
    ignore_column::Integer=-1
) where {I, C <: Union{CoeffZp, CompositeCoeffZp}, A <: Union{CoeffZp, CompositeCoeffZp}}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    n_nonzeros = 0
    new_pivot_column = -1

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
        if !isassigned(pivots, i) || ignore_column == i
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        linalg_vector_addmul_sparsedense!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    if n_nonzeros == 0
        return true
    end

    if end_column != length(row)
        linalg_dense_row_mod_p!(row, arithmetic, end_column + 1, length(row))
    end

    # form the resulting row in sparse format
    resize!(new_sparse_row_support, n_nonzeros)
    resize!(new_sparse_row_coeffs, n_nonzeros)
    linalg_extract_sparse_row!(
        new_sparse_row_support,
        new_sparse_row_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

# Generic fallback.
function linalg_reduce_dense_row_by_pivots_sparse!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmetic,
    active_reducers=nothing;
    ignore_column::Integer=-1
) where {I, C <: Coeff}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)

    n_nonzeros = 0
    new_pivot_column = -1

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || ignore_column == i
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        # Locate the support and the coefficients of the reducer row
        indices = pivots[i]
        if i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(indices) == length(coeffs)

        linalg_vector_addmul_sparsedense!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if n_nonzeros == 0
        return true
    end

    # form the resulting row in sparse format
    resize!(new_sparse_row_support, n_nonzeros)
    resize!(new_sparse_row_coeffs, n_nonzeros)
    linalg_extract_sparse_row!(
        new_sparse_row_support,
        new_sparse_row_coeffs,
        row,
        convert(Int, start_column),
        ncols
    )

    false
end

###
# Basic routines for dense and sparse vectors

function linalg_new_empty_sparse_row(::Type{C}) where {C <: Coeff}
    Vector{ColumnLabel}(), Vector{C}()
end

function linalg_dense_row_mod_p!(
    row::Vector{T},
    arithmetic::AbstractArithmeticZp{T},
    from::Int=1,
    to::Int=length(row)
) where {T <: Union{CoeffZp, CompositeCoeffZp}}
    # This loop is usually successfully vectorized for all kinds of arithmetic.
    # Still, the resulting code may be not optimal
    @fastmath @inbounds for i in from:to
        row[i] = mod_p(row[i], arithmetic)
    end
    nothing
end

function linalg_row_make_monic!(
    row::Vector{T},
    arithmetic::AbstractArithmeticZp{A, T},
    first_nnz_index::Int=1
) where {A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant !iszero(row[first_nnz_index])

    @inbounds lead = row[first_nnz_index]
    isone(lead) && return lead

    @inbounds pinv = inv_mod_p(A(lead), arithmetic) % T
    @inbounds row[first_nnz_index] = one(T)
    @inbounds for i in (first_nnz_index + 1):length(row)
        row[i] = mod_p(A(row[i]) * A(pinv), arithmetic) % T
    end

    pinv
end

# Generic fallback
function linalg_row_make_monic!(
    row::Vector{T},
    arithmetic::AbstractArithmetic,
    first_nnz_index::Int=1
) where {T <: Coeff}
    @invariant !iszero(row[first_nnz_index])

    @inbounds lead = row[first_nnz_index]
    isone(lead) && return lead

    @inbounds pinv = inv(lead)
    @inbounds row[1] = one(row[1])
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end

    pinv
end

# Linear combination of dense vector and sparse vector
# with explicit reduction modulo a prime.
function linalg_vector_addmul_sparsedense_mod_p!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticZp
) where {I, A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

    @inbounds mul = divisor(arithmetic) - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        # if A === T, then the type cast is a no-op
        a = row[idx] + A(mul) * A(coeffs[j])
        row[idx] = mod_p(a, arithmetic)
    end

    mul
end

# Linear combination of dense vector and sparse vector.
function linalg_vector_addmul_sparsedense!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticZp
) where {I, A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

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
function linalg_vector_addmul_sparsedense!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::SignedArithmeticZp{A, T}
) where {I, A <: CoeffZp, T <: CoeffZp}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

    p2 = arithmetic.p2
    @invariant row[indices[1]] < typemax(T)
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
function linalg_vector_addmul_sparsedense!(
    row::Vector{CompositeNumber{N, A}},
    indices::Vector{I},
    coeffs::Vector{CompositeNumber{N, T}},
    arithmetic::SignedCompositeArithmeticZp{CompositeNumber{N, A}, CompositeNumber{N, T}}
) where {I, A <: CoeffZp, T <: CoeffZp, N}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

    p2 = arithmetic.p2s
    @invariant row[indices[1]] < typemax(T)
    @inbounds mul = row[indices[1]].data .% T

    @fastmath @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx].data .- A.(mul) .* A.(coeffs[j].data)
        a = a .+ ((a .>> (8 * sizeof(A) - 1)) .& p2.data)
        row[idx] = CompositeNumber(a)
    end

    nothing
end

# Generic fallback.
function linalg_vector_addmul_sparsedense!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmetic
) where {I, A <: Coeff, T <: Coeff}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + A(mul) * coeffs[j]
    end

    nothing
end

# Loads the given sparse row into the given dense row.
# Zero out the entries of the dense row before that.
function linalg_load_sparse_row!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T}
) where {I, A <: Coeff, T <: Coeff}
    @invariant length(indices) == length(coeffs)

    @inbounds z = zero(row[1])
    @inbounds for i in 1:length(row)
        row[i] = z
    end

    @inbounds for j in 1:length(indices)
        row[indices[j]] = A(coeffs[j])
    end

    nothing
end

# Extracts nonzeroes from the given dense row into the given sparse row.
# Returns the number of extracted nonzeroes.
function linalg_extract_sparse_row!(
    indices::Vector{I},
    coeffs::Vector{T},
    row::Vector{A},
    from::Int,
    to::Int
) where {I, T <: Coeff, A <: Coeff}
    # NOTE: assumes that the sparse row has the necessary capacity
    @invariant length(indices) == length(coeffs)
    @invariant 1 <= from <= to <= length(row)

    nnz = 1
    @inbounds for i in from:to
        if !iszero(row[i])
            indices[nnz] = i
            @invariant row[i] <= typemax(T)
            coeffs[nnz] = row[i] % T
            nnz += 1
        end
    end

    @invariant nnz - 1 <= length(coeffs)
    nnz - 1
end

# Generic fallback.
function linalg_extract_sparse_row!(
    indices::Vector{I},
    coeffs::Vector{T},
    row::Vector{T},
    from::Int,
    to::Int
) where {I, T <: Coeff}
    # NOTE: assumes that the sparse row has the necessary capacity
    @invariant length(indices) == length(coeffs)
    @invariant 1 <= from <= to <= length(row)

    nnz = 1
    @inbounds for i in from:to
        if !iszero(row[i])
            indices[nnz] = i
            coeffs[nnz] = row[i]
            nnz += 1
        end
    end

    @invariant nnz - 1 <= length(coeffs)
    nnz - 1
end
