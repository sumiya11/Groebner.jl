# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_deterministic_sparse_threaded!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "linalg_deterministic_sparse!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    linalg_reduce_matrix_lower_part_threaded_cas!(matrix, basis, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)
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
function linalg_reduce_matrix_lower_part_threaded_cas!(
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
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    resize!(matrix.sentinels, ncols)

    # If there is a pivot with the leading column i, then we set
    #   sentinels[i] = 1.
    # If sentinels[i] = 1, and the pivot is synced, then we set
    #   sentinels[i] = 2.
    # Otherwise, sentinels[i] = 0.
    # Once sentinels[i] becomes 1 or 2, this equality is unchanged.
    sentinels = matrix.sentinels
    @inbounds for i in 1:ncols
        sentinels[i] = 0
    end
    @inbounds for i in 1:nup
        sentinels[matrix.upper_rows[i][1]] = 1
    end

    # Allocate thread-local buffers
    buffers_row = map(_ -> zeros(AccumType, ncols), 1:nthreads())

    # NOTE: by default, @threads uses the :dynamic execution schedule, which
    # does not guarantee that threadid() is constant within one iteration
    @inbounds Base.Threads.@threads for i in 1:nlow
        t_id = threadid()
        t_local_row = buffers_row[t_id]
        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)

        # Select a row to be reduced from the lower part of the matrix
        # NOTE: no copy of coefficients is needed
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the coefficients into a dense array
        linalg_load_sparse_row!(t_local_row, sparse_row_support, sparse_row_coeffs)

        first_nnz_column = sparse_row_support[1]

        # In a somewhat ideal scenario, the while loop below executes just once
        success = false
        @inbounds while !success
            @invariant 1 <= first_nnz_column <= length(t_local_row)

            # Note that this call may only mutate the first three arguments, and
            # the arguments do not overlap in memory with each other
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
                tmp_pos=-1
            )

            # If the row is fully reduced
            zeroed && break

            # At this point, we have discovered a row that is not reduced to
            # zero (yet). Thus, two outcomes are possible: 
            # 1. sentinels[...] is 0, meaning that this row is potentially a new
            #   pivot row. We try to set sentinels[...] to 1, update the pivots
            #   and the matrix, and break out of the loop.
            #   If assigning sentinels[...] to 1 fails, go to 2.
            # 2. sentinels[...] is already 1, meaning that this row can be
            #   further reduced on the following iteration

            # Sync point. Everything before this point becomes visible to
            # all other threads once they reach this point.
            # NOTE: Note the absense of a total ordering on atomic operations
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

                linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)

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
                #
                # Roughly speaking, in single-thread, this construction is
                # similar to an explicit barrier in C
                #   asm volatile("" ::: "memory");
                # which can be used to tell the compiler that the writes and
                # loads cannot be freely moved below or above the barrier.
                Atomix.set!(
                    Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                    Int8(1),
                    Atomix.release
                )
            else
                # go for another iteration
                first_nnz_column = new_sparse_row_support[1]
            end
        end

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    true
end

function linalg_reduce_dense_row_by_pivots_sparse_threadsafe0!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmeticZp{A, C},
    sentinels::Vector{Int8},
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

        # Say, there are two threads, A and B. Thread A non-atomically writes to
        # pivots[i]. Then loading pivots[i] in thread B may produce a value out
        # of thin air.
        #
        # Modifying pivots[i] atomically is impossible in Julia.
        # Thus, we may need additional syncronization in this case.
        # (Note that this is usually not an issue for x86-64, since it performs
        # well-aligned writes and loads atomically by default, but we should not
        # rely on that)
        # 
        # From https://en.cppreference.com/w/cpp/atomic/memory_order, if thread
        # B performs (at least) acquire atomic memory load, then "all memory
        # writes that happened-before the atomic store from the point of view of
        # thread A, become visible side-effects in thread B".
        #
        # In other words, consider the following "sequence" of writes/loads:
        # thread A:
        #   1. write pivots[i] = X
        #   2. write sentinels[i] = 1 (atomic)
        # thread B:
        #   3. load sentinels[i] (atomic)
        #   4. load pivots[i]
        # If 2. happened before 3., then 1. happened before 4.
        #
        # For this reason, we perform the atomic write 2. in
        # `linalg_reduce_matrix_lower_part_threaded_cas_maybe_correct!`, and
        # perform the atomic load 3. here.
        if iszero(Atomix.get(Atomix.IndexableRef(sentinels, (i,)), Atomix.acquire))
            if new_pivot_column == -1
                new_pivot_column = i
            end
            n_nonzeros += 1
            continue
        end
        # At this point, the value of pivots[i] is in sync with any other thread
        # that modifies, modified, or will modify it

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

@timeit function linalg_reduce_matrix_lower_part_threaded_task_cas!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Calculate the size of the block
    nblocks = nthreads()
    nblocks = min(nblocks, nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    @log level = -3 "" nblocks rem rowsperblock nlow

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
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

function _reduce_matrix_lower_part_threaded_task_cas_worker!(
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
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    sentinels = matrix.sentinels

    @inbounds for i in block_start:block_end
        # Select the row from the lower part of the matrix to be reduced
        sparse_row_support = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the coefficients into a dense array
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        success = false
        while !success

            # Reduce the row with respect to the known `pivots` from the upper part
            # of the matrix.
            # NOTE: this also does partial interreduction of the lower matrix rows.
            # Upon discovery of a new pivot from the lower part of the matrix, we
            # add the pivot to the pool of available pivots
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
            zeroed && break

            @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)

            # Set a reference to the coefficients of this row in the matrix
            old, success = Atomix.replace!(
                Atomix.IndexableRef(sentinels, (Int(new_sparse_row_support[1]),)),
                Int8(0),
                Int8(1),
                Atomix.acquire_release,
                Atomix.acquire
            )

            if success
                @invariant iszero(old)
                linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)
                matrix.some_coeffs[i] = new_sparse_row_coeffs
                matrix.lower_to_coeffs[new_sparse_row_support[1]] = i
                pivots[new_sparse_row_support[1]] = new_sparse_row_support
            else
                sparse_row_support = new_sparse_row_support
                sparse_row_coeffs = new_sparse_row_coeffs
            end
        end

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end
end

@timeit function linalg_reduce_matrix_lower_part_threaded_task_lock_free!(
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
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
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
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        if iszero(matrix.sentinels[i])
            continue
        end

        # Select the row from the lower part of the matrix to be reduced
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = matrix.some_coeffs[i]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the coefficients into a dense array
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add the pivot to the pool of available pivots
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

        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)
        linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)

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

function _reduce_matrix_lower_part_threaded_task_lock_free_worker!(
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
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in block_start:block_end
        # Select the row from the lower part of the matrix to be reduced
        sparse_row_support = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Load the coefficients into a dense array
        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add the pivot to the pool of available pivots
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

        @invariant length(new_sparse_row_coeffs) == length(new_sparse_row_support)
        linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)

        # Store the new row in the matrix, AND add it to the active pivots
        matrix.some_coeffs[i] = new_sparse_row_coeffs
        matrix.lower_rows[i] = new_sparse_row_support
        matrix.sentinels[i] = 1

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end
end
