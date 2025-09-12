# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_randomized_sparse_threaded!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    tasks::Int,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    # Reduce CD with AB
    linalg_randomized_reduce_matrix_lower_part_threaded_cas!(matrix, basis, tasks, arithmetic, rng)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

###
# Low level

function _linalg_randomized_reduce_matrix_lower_part_threaded_cas!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    chunk::Vector{Int},
    pivots::Vector{Vector{ColumnLabel}},
    row_index_to_coeffs::Vector{Int},
    sentinels::Vector{Int8},
    rowsperblock::Int,
    arithmetic::AbstractArithmetic,
    t_local_row::Vector{AccumType},
    t_local_rows_multipliers::Vector{AccumType},
    t_local_rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, nlow = matrix_nrows_filled(matrix)
    _, ncols = size(matrix)
    for i in chunk
        nrowsupper = min(i * rowsperblock, nlow)
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

        new_pivots_count = 0
        @inbounds while new_pivots_count < nrowstotal
            for j in 1:nrowstotal
                t_local_rows_multipliers[j] = mod_p(rand(t_local_rng, AccumType), arithmetic)
            end
            t_local_row .= AccumType(0)
            first_nnz_col = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                sparse_row_support = matrix.lower_rows[rowidx]
                sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[rowidx]]
                @invariant length(sparse_row_support) == length(sparse_row_coeffs)

                first_nnz_col = min(first_nnz_col, sparse_row_support[1])
                for l in 1:length(sparse_row_coeffs)
                    colidx = sparse_row_support[l]
                    t_local_row[colidx] = mod_p(
                        t_local_row[colidx] +
                        AccumType(t_local_rows_multipliers[k]) * AccumType(sparse_row_coeffs[l]),
                        arithmetic
                    )
                end
            end

            # In a somewhat ideal scenario, the loop below executes just once
            success = false
            while !success
                # Reduce the combination by pivots
                @invariant 1 <= first_nnz_col <= length(t_local_row)
                zeroed = linalg_reduce_dense_row_by_pivots_sparse_threadsafe0!(
                    new_sparse_row_support,
                    new_sparse_row_coeffs,
                    t_local_row,
                    matrix,
                    basis,
                    pivots,
                    first_nnz_col,
                    ncols,
                    arithmetic,
                    sentinels
                )

                if zeroed
                    # set the flag so that the outer loop exits immediately
                    new_pivots_count = nrowstotal
                    break
                end

                # Sync point. Everything before this point becomes visible to
                # other threads once they reach this point.
                # Note the absense of a total order on atomic operations
                old, success = Atomix.replace!(
                    Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                    Int8(0),
                    Int8(1),
                    Atomix.acquire_release,
                    Atomix.acquire
                )

                @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)

                if success
                    # new pivot is formed, wake up from the loop
                    @invariant iszero(old)

                    new_pivots_count += 1
                    absolute_row_index = (i - 1) * rowsperblock + new_pivots_count

                    linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

                    matrix.some_coeffs[absolute_row_index] = new_sparse_row_coeffs
                    matrix.lower_to_coeffs[new_sparse_row_support[1]] = absolute_row_index
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
                    
                    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
                else
                    # go for another iteration
                    first_nnz_col = new_sparse_row_support[1]
                end
            end
        end
    end
end

function linalg_randomized_reduce_matrix_lower_part_threaded_cas!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    tasks::Int,
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    nblocks = linalg_nblocks_in_randomized(nlow)
    nblocks = max(nblocks, min(nlow, tasks))
    @assert nblocks <= nlow
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
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
    buffers_rows_multipliers = map(_ -> zeros(AccumType, rowsperblock), 1:tasks)
    buffers_row = map(_ -> zeros(AccumType, ncols), 1:tasks)
    buffers_rng = map(_ -> copy(rng), 1:tasks)

    tasks = min(tasks, nblocks)
    data_chunks = split_round_robin(1:nblocks, tasks)
    task_results = Vector{Task}(undef, tasks)
    for (tid, chunk) in enumerate(data_chunks)
        t_local_rows_multipliers = buffers_rows_multipliers[tid]
        t_local_row = buffers_row[tid]
        t_local_rng = buffers_rng[tid]
        task = @spawn begin
            _linalg_randomized_reduce_matrix_lower_part_threaded_cas!(
                matrix,
                basis,
                chunk,
                pivots,
                row_index_to_coeffs,
                sentinels,
                rowsperblock,
                arithmetic,
                t_local_row,
                t_local_rows_multipliers,
                t_local_rng
            )
        end
        task_results[tid] = task
    end
    for task in task_results
        wait(task)
    end

    true
end
