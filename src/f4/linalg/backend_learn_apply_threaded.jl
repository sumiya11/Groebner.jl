# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_learn_sparse_threaded!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    tasks::Int,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    # Reduce CD with AB
    linalg_learn_reduce_matrix_lower_part_threaded!(trace, matrix, basis, tasks, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function _linalg_learn_reduce_matrix_lower_part_threaded!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    chunk::Vector{Int},
    pivots::Vector{Vector{ColumnLabel}},
    row_idx_to_coeffs::Vector{Int},
    sentinels::Vector{Int8},
    arithmetic::AbstractArithmetic,
    t_local_row::Vector{AccumType},
    t_local_not_reduced_to_zero::Vector{Int},
    t_local_useful_reducers::Set{Int},
    t_local_reducer_rows::Vector{Int}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    @inbounds for i in chunk
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        linalg_load_sparse_row!(t_local_row, sparse_row_support, sparse_row_coeffs)

        # Additionally record the indices of rows that participated in reduction
        # of the given row
        empty!(t_local_reducer_rows)
        first_nnz_column = sparse_row_support[1]

        success = false
        @inbounds while !success
            @invariant 1 <= first_nnz_column <= length(t_local_row)

            zeroed = linalg_reduce_dense_row_by_pivots_sparse_threadsafe0!(
                new_sparse_row_support,
                new_sparse_row_coeffs,
                t_local_row,
                matrix,
                basis,
                pivots,
                first_nnz_column,
                ncols,
                arithmetic,
                sentinels,
                t_local_reducer_rows
            )

            # if fully reduced
            zeroed && break

            # Sync point. Everything before this point becomes visible to
            # all other threads once they reach this point.
            # Note the absense of a total ordering on atomic operations
            old, success = Atomix.replace!(
                Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                Int8(0),
                Int8(1),
                Atomix.acquire_release,
                Atomix.acquire
            )

            @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)

            if success
                # new pivot is formed, wake up from the inner loop
                @invariant iszero(old)

                # NOTE: we are not recording reducers from the lower part of the matrix
                push!(t_local_not_reduced_to_zero, i)
                for reducer_row in t_local_reducer_rows
                    push!(t_local_useful_reducers, reducer_row)
                end

                linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

                matrix.some_coeffs[i] = new_sparse_row_coeffs
                matrix.lower_to_coeffs[new_sparse_row_support[1]] = i
                pivots[new_sparse_row_support[1]] = new_sparse_row_support

                # This atomic write is paired by an atomic load in
                # `linalg_reduce_dense_row_by_pivots_sparse_threadsafe0!`. The
                # purpose of these atomics is to fence the above assignment
                #    pivots[new_sparse_row_support[1]] = new_sparse_row_support
                # 
                # This way, the next atomic load from the same location in
                # `sentinels` will also receive all relevant modifications to
                # `pivots` as a side-effect.
                Atomix.set!(
                    Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                    Int8(2),
                    Atomix.release
                )

                new_sparse_row_support, new_sparse_row_coeffs =
                    linalg_new_empty_sparse_row(CoeffType)
            else
                # go for another iteration
                first_nnz_column = new_sparse_row_support[1]
            end
        end
    end
end

function linalg_learn_reduce_matrix_lower_part_threaded!(
    trace::Trace,
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    tasks::Int,
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.sentinels, ncols)

    # If there is a pivot with the leading column i, then we set
    #   sentinels[i] = 1.
    # If sentinels[i] = 1, and if the pivot is synced, then we set
    #   sentinels[i] = 2.
    # Otherwise, sentinels[i] = 0.
    # Once sentinels[i] becomes 1 or 2, this is unchanged.
    sentinels = matrix.sentinels
    @inbounds for i in 1:ncols
        sentinels[i] = 0
    end
    @inbounds for i in 1:nup
        sentinels[matrix.upper_rows[i][1]] = 2
    end

    # Allocate thread-local buffers
    buffers_row = map(_ -> zeros(AccumType, ncols), 1:tasks)
    not_reduced_to_zero = map(_ -> Vector{Int}(), 1:tasks)
    useful_reducers = map(_ -> Set{Int}(), 1:tasks)
    reducer_rows = map(_ -> Vector{Int}(), 1:tasks)

    tasks = min(tasks, nlow)
    data_chunks = split_round_robin(1:nlow, tasks)
    task_results = Vector{Task}(undef, tasks)
    for (tid, chunk) in enumerate(data_chunks)
        t_local_row = buffers_row[tid]
        t_local_not_reduced_to_zero = not_reduced_to_zero[tid]
        t_local_useful_reducers = useful_reducers[tid]
        t_local_reducer_rows = reducer_rows[tid]
        task = @spawn begin
            _linalg_learn_reduce_matrix_lower_part_threaded!(
                matrix,
                basis,
                chunk,
                pivots,
                row_idx_to_coeffs,
                sentinels,
                arithmetic,
                t_local_row,
                t_local_not_reduced_to_zero,
                t_local_useful_reducers,
                t_local_reducer_rows
            )
        end
        task_results[tid] = task
    end
    for task in task_results
        wait(task)
    end

    # Update the tracer information
    # NOTE: we sort reducers by their original position in the array of pivots.
    # This way, the rows are already sorted at the apply stage.
    useful_reducers_sorted = sort(collect(reduce(union, useful_reducers)))
    push!(
        trace.matrix_infos,
        (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
    )
    _not_reduced_to_zero = sort(reduce(vcat, not_reduced_to_zero))
    push!(trace.matrix_nonzeroed_rows, _not_reduced_to_zero)
    push!(
        trace.matrix_upper_rows,
        (
            map(f -> matrix.upper_to_coeffs[f], useful_reducers_sorted),
            map(f -> matrix.upper_to_mult[f], useful_reducers_sorted)
        )
    )
    push!(
        trace.matrix_lower_rows,
        (row_idx_to_coeffs[_not_reduced_to_zero], matrix.lower_to_mult[_not_reduced_to_zero])
    )

    true
end
