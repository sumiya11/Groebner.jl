# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_learn_sparse_threaded!(
        trace::TraceF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_learn_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Reduce CD with AB
    linalg_learn_reduce_matrix_lower_part_threaded!(trace, matrix, basis, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

###
# Low level

function linalg_learn_reduce_matrix_lower_part_threaded!(
        trace::TraceF4,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    if nthreads() == 1
        @log :info """
        Using multi-threaded linear algebra with nthreads() == 1. 
        Something probably went wrong."""
    end
    @assert false # this may be broken

    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
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
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        t_id = threadid()
        row = row_buffers[t_id]
        t_alloc = ctx.t_alloc[t_id]
        useful_reducers = useful_reducers_buffers[t_id]
        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        success = false
        while !success
            first_nnz_col = sparse_row_support[1]

            # Reduce the combination by rows from the upper part of the matrix
            n_nonzeros = linalg_reduce_dense_row_by_pivots_sparse!(
                new_sparse_row_support,
                new_sparse_row_coeffs,
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

            @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)

            old, success = Atomix.replace!(
                Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                Int8(0),
                Int8(1),
                Atomix.acquire_release,
                Atomix.acquire
            )

            if success
                @invariant iszero(old)

                linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)
                matrix.some_coeffs[i] = new_sparse_row_coeffs
                matrix.lower_to_coeffs[new_sparse_row_support[1]] = i
                pivots[new_sparse_row_support[1]] = new_sparse_row_support

                # NOTE: we are not recording reducers from the lower part of the matrix
                not_reduced_to_zero[i] = 1
            else
                sparse_row_support = new_sparse_row_support
                sparse_row_coeffs = new_sparse_row_coeffs
            end
        end

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
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
