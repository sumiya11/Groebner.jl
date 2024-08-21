# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_learn_sparse!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_learn_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Reduce CD with AB
    linalg_learn_reduce_matrix_lower_part!(trace, matrix, basis, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function linalg_apply_sparse!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    # NOTE: here, we do not need to sort the rows in the upper part, as they
    # have already been collected in the right order
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_apply_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Reduce CD with AB
    flag = linalg_apply_reduce_matrix_lower_part!(trace, matrix, basis, arithmetic)
    if !flag
        return flag
    end
    # Interreduce CD
    linalg_apply_interreduce_matrix_pivots!(trace, matrix, basis, arithmetic)

    true
end

function linalg_learn_deterministic_sparse_interreduction!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log :matrix "linalg_learn_deterministic_sparse_interreduction!"
    @log :matrix matrix_string_repr(matrix)

    # Prepare the matrix
    linalg_prepare_matrix_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    linalg_learn_interreduce_matrix_pivots!(
        trace,
        matrix,
        basis,
        arithmetic,
        reversed_rows=true
    )

    true
end

function linalg_apply_deterministic_sparse_interreduction!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    @log :matrix "linalg_apply_deterministic_sparse_interreduction!"
    @log :matrix matrix_string_repr(matrix)

    # Prepare the matrix
    linalg_prepare_matrix_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    flag = linalg_apply_interreduce_matrix_pivots!(
        trace,
        matrix,
        basis,
        arithmetic,
        reversed_rows=true
    )

    flag
end

###
# Low level

record_active_reducer(active_reducers::Nothing, matrix, idx) = nothing
function record_active_reducer(active_reducers, matrix, idx)
    push!(active_reducers, (idx, matrix.upper_to_coeffs[idx], matrix.upper_to_mult[idx]))
    nothing
end

# Returns `false` if any row reduced to zero (since we expect that on the apply
# stage the rows are linearly independent)
function linalg_apply_reduce_matrix_lower_part!(
    trace::Trace,
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
    row = zeros(AccumType, ncols)
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        # Select the row from the lower part of the matrix to be reduced
        sparse_row_support = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]

        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load coefficients into a dense array
        # TODO!!!: if `row` was fully reduced to zero on the previous iteration,
        # then do not set it to zero in here
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
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
        if zeroed
            return false
        end

        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)
        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    true
end

function linalg_learn_interreduce_matrix_pivots!(
    trace::Trace,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    # Perform interreduction
    flag, _, not_reduced_to_zero = linalg_interreduce_matrix_pivots!(
        matrix,
        basis,
        arithmetic,
        reversed_rows=reversed_rows
    )
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
    # TODO: see "TODO: (I)" in src/groebner/groebner.jl
    push!(trace.matrix_lower_rows, (Vector{Int}(), Vector{Int}()))

    true
end

function linalg_apply_interreduce_matrix_pivots!(
    trace::Trace,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    flag, any_zeroed, _ = linalg_interreduce_matrix_pivots!(
        matrix,
        basis,
        arithmetic,
        reversed_rows=reversed_rows
    )
    flag && !any_zeroed
end

function linalg_learn_reduce_matrix_lower_part!(
    trace::Trace,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    not_reduced_to_zero = Vector{Int}()
    pivot_indices = Vector{Int}()
    useful_reducers = Set{Tuple{Int, Int, MonomId}}()
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    reducer_rows = Tuple{Int, Int, MonomId}[]

    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Additionally record the indices of rows that participated in reduction
        # of the given row
        empty!(reducer_rows)
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
            reducer_rows,
            tmp_pos=-1
        )

        # if fully reduced
        zeroed && continue

        # NOTE: we are not recording reducers from the lower part of the matrix
        push!(not_reduced_to_zero, i)
        push!(pivot_indices, new_sparse_row_support[1])
        for reducer_row in reducer_rows
            push!(useful_reducers, reducer_row)
        end

        @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)
        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
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
    # push!(trace.matrix_pivot_indices, map(sgn -> matrix.column_to_monom[sgn], pivot_indices))

    true
end
