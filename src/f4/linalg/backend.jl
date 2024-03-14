# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

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

    @log :matrix "linalg_deterministic_sparse!"
    @log :matrix matrix_string_repr(matrix)

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
    @log :matrix "linalg_deterministic_sparse_interreduction!"
    @log :matrix matrix_string_repr(matrix)

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
@timeit function linalg_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    # TODO: allocate just once, and then reuse for later iterations of F4
    row = zeros(AccumType, ncols)
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
            arithmetic,
            tmp_pos=-1
        )

        # If the row is fully reduced
        zeroed && continue

        @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)
        linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)

        # Store the new row in the matrix, AND add the new row to the list of
        # known pivots
        pivots[new_sparse_row_support[1]] = new_sparse_row_support

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
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
@timeit function linalg_interreduce_matrix_pivots!(
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

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
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
            matrix.lower_rows[new_pivots] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] =
                new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[new_pivots]
        else
            matrix.lower_rows[nupper - new_pivots + 1] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] =
                new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[nupper - new_pivots + 1]
        end
    end

    matrix.npivots = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    true, any_zeroed, not_reduced_to_zero
end

# Puts the AB part of the matrix in the RREF, inplace.
@timeit function linalg_interreduce_matrix_upper_part!(
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
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in nup:-1:1
        # Locate the support and the coefficients of the row from the upper part
        # of the matrix
        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]

        # Extract the coefficients into a dense array
        @invariant isone(sparse_row_coeffs[1])
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: note the `tmp_pos=first_nnz_column` argument. It ensures that we
        # do not reduce the row with itself. Maybe think of a better name?..
        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
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

        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)
        linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)

        # Update the support and the coefficients of the vector
        matrix.upper_coeffs[i] = new_sparse_row_coeffs
        matrix.upper_rows[i] = new_sparse_row_support

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    matrix.upper_part_is_rref = true

    true
end

function linalg_reduce_matrix_lower_part_invariant_pivots!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
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
            arithmetic,
            tmp_pos=-1
        )

        # NOTE: no normalization here
        # NOTE: the array may be empty, signalling the row is zero
        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        matrix.lower_rows[i] = new_sparse_row_support
        matrix.lower_to_coeffs[i] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    matrix.npivots = matrix.nrows_filled_lower = matrix.nrows_filled_lower
end

function linalg_reduce_matrix_lower_part_any_nonzero!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
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

function linalg_prepare_matrix_pivots_in_interreduction!(
    matrix::MacaulayMatrix,
    basis::Basis
)
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

# The most generic specialization.
# 
# Reduces the given dense `row` by the given `pivots`.
# It is assumed that the leading coefficients of pivots are equal to one.
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
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
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
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # At this point, we have determined that pivots[i] is a valid pivot.

        # Locate the support and the coefficients of the pivot row
        pivot_support = pivots[i]
        if exact_column_mapping
            # if pivot comes from the new pivots
            pivot_coeffs = matrix.some_coeffs[tmp_pos]
        elseif i <= nleft
            # if pivot comes from the original pivots
            if matrix.upper_part_is_rref || computing_rref
                pivot_coeffs = matrix.upper_coeffs[i]
            else
                pivot_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            end
            record_active_reducer(active_reducers, matrix, i)
        else
            pivot_coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(pivot_support) == length(pivot_coeffs)

        linalg_vector_addmul_sparsedense_mod_p!(
            row,
            pivot_support,
            pivot_coeffs,
            arithmetic
        )

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
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
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
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
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
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false,
    computing_rref::Bool=false
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
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
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

        linalg_vector_addmul_sparsedense!(row, indices, coeffs, arithmetic)

        @invariant iszero(row[i])
    end

    # all reduced to zero!
    if n_nonzeros == 0
        return true
    end

    if end_column != length(row)
        linalg_dense_row_mod_p!(row, arithmetic, end_columns + 1, length(row))
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

# Specialization for arithmetic in the rational numbers
function linalg_reduce_dense_row_by_pivots_sparse!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
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

    n_nonzeros = 0
    new_pivot_column = -1

    @inbounds for i in start_column:end_column
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
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
    # NOTE: This loop is usually successfully vectorized by the compiler for all
    # kinds of arithmetic. Still, the resulting native code may be not optimal
    @fastmath @inbounds for i in from:to
        row[i] = mod_p(row[i], arithmetic)
    end
    nothing
end

# Normalize the row to have the leading coefficient equal to 1
function linalg_normalize_row!(
    row::Vector{T},
    arithmetic::AbstractArithmeticZp{A, T},
    first_nnz_index::Int=1
) where {A <: Union{CoeffZp, CompositeCoeffZp}, T <: Union{CoeffZp, CompositeCoeffZp}}
    @invariant !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    isone(lead) && return lead

    @inbounds pinv = inv_mod_p(A(lead), arithmetic) % T
    @inbounds row[first_nnz_index] = one(T)
    @inbounds for i in (first_nnz_index + 1):length(row)
        row[i] = mod_p(A(row[i]) * A(pinv), arithmetic) % T
    end

    pinv
end

# Normalize the row to have the leading coefficient equal to 1
function linalg_normalize_row!(
    row::Vector{T},
    arithmetic::AbstractArithmeticQQ{T},
    first_nnz_index::Int=1
) where {T <: CoeffQQ}
    @invariant !iszero(row[first_nnz_index])

    lead = row[first_nnz_index]
    isone(lead) && return lead

    @inbounds pinv = inv(lead)
    @inbounds row[1] = one(T)
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end

    pinv
end

# Linear combination of dense vector and sparse vector modulo a prime.
# The most generic version.
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
# The most generic version.
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
# Specialization for SignedArithmeticZp.
function linalg_vector_addmul_sparsedense!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::SignedArithmeticZp{A, T}
) where {I, A <: CoeffZp, T <: CoeffZp}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

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
function linalg_vector_addmul_sparsedense!(
    row::Vector{CompositeInt{N, A}},
    indices::Vector{I},
    coeffs::Vector{CompositeInt{N, T}},
    arithmetic::SignedCompositeArithmeticZp{CompositeInt{N, A}, CompositeInt{N, T}}
) where {I, A <: CoeffZp, T <: CoeffZp, N}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

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
function linalg_vector_addmul_sparsedense!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::AbstractArithmeticQQ{T}
) where {I, T <: CoeffQQ}
    @invariant isone(coeffs[1])
    @invariant length(indices) == length(coeffs)
    @invariant !isempty(indices)

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
    end

    nothing
end

# Linear combination of two dense vectors modulo a prime
function linalg_vector_addmul_densedense_mod_p!(
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
function linalg_vector_addmul_densedense!(
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
function linalg_load_sparse_row!(
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
function linalg_load_sparse_row!(
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
function linalg_extract_sparse_row!(
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
function linalg_extract_sparse_row!(
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
