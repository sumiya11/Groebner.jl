# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

# These are experimental linear algebra backends. 
# These are not used and not tested.

###
# High level

function linalg_randomized_hashcolumns_sparse!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_randomized_hashcolumns_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Reduce CD with AB
    linalg_randomized_hashcolumns_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

function linalg_direct_rref_sparse!(
        matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_direct_rref_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Produce the RREF of AB
    linalg_interreduce_matrix_upper_part!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    linalg_randomized_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

function linalg_direct_rref_sparsedense!(
        matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_direct_rref_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Produce the RREF of AB
    linalg_interreduce_matrix_upper_part_sparsedense!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    linalg_reduce_matrix_lower_part_sparsedense!(matrix, basis, arithmetic)
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)

    true
end

###
# Low level

@timeit function linalg_reduce_matrix_lower_part_sparsedense!(
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
        sparse_row_support = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        sparse_row_coeffs = basis.coeffs[matrix.lower_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        # Extract the sparse row into two dense arrays, so that the row equals
        # [row1..., row2...]
        for j in 1:nleft
            row1[j] = zero(C)
        end
        for j in 1:nright
            row2[j] = zero(C)
        end
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                row1[idx] = sparse_row_coeffs[j]
            else
                row2[idx - nleft] = sparse_row_coeffs[j]
            end
        end

        # Reduce the rows simultaneously
        first_nnz_column = sparse_row_support[1]
        pivot = linalg_reduce_dense_row_by_pivots_sparsedense!(
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

        linalg_row_make_monic!(row2, arithmetic, pivot)
        matrix.D_coeffs_dense[i] = row2
        matrix.pivot_indices[i] = pivot
    end

    true
end

@timeit function linalg_interreduce_matrix_upper_part_sparsedense!(
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

        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]

        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        @invariant isone(sparse_row_coeffs[1])

        # Extract the sparse row into two dense arrays, so that the row equals
        # [row1..., row2...]
        for j in 1:nleft
            row1[j] = zero(AccumType)
        end
        for j in 1:nright
            row2[j] = zero(AccumType)
        end
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                row1[idx] = sparse_row_coeffs[j]
            else
                row2[idx - nleft] = sparse_row_coeffs[j]
            end
        end

        first_nnz_column = sparse_row_support[1]
        _ = linalg_reduce_dense_row_by_pivots_sparsedense!(
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
        linalg_row_make_monic!(row2, arithmetic, first_nnz_column)

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
            linalg_vector_addmul_densedense!(to_be_reduced, row2, mul, arithmetic)
        end

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
        resize!(new_sparse_row_support, length(row2))
        resize!(new_sparse_row_coeffs, length(row2))
        cnt = linalg_extract_sparse_row!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row2,
            first_nnz_column,
            nright
        )
        resize!(new_sparse_row_support, cnt)
        for j in 1:length(new_sparse_row_support)
            new_sparse_row_support[j] = nleft + new_sparse_row_support[j]
        end
        resize!(new_sparse_row_coeffs, cnt)

        first_nnz_column = first_nnz_column + nleft
        @invariant isone(new_sparse_row_coeffs[1])

        new_pivots += 1
        not_reduced_to_zero[new_pivots] = first_nnz_column
        # update row and coeffs
        if !reversed_rows
            matrix.lower_rows[new_pivots] = new_sparse_row_support
            matrix.lower_to_coeffs[first_nnz_column] = new_pivots
            matrix.some_coeffs[new_pivots] = new_sparse_row_coeffs
        else
            matrix.lower_rows[nlower - new_pivots + 1] = new_sparse_row_support
            matrix.lower_to_coeffs[first_nnz_column] = new_pivots
            matrix.some_coeffs[new_pivots] = new_sparse_row_coeffs
        end
    end

    matrix.npivots = matrix.nrows_filled_lower = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    true, any_zeroed, not_reduced_to_zero
end

function linalg_randomized_hashcolumns_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)
    nleft, nright = matrix_ncols_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    indir_up_left, indir_up_right =
        Vector{Vector{Int}}(undef, nup), Vector{Vector{Int}}(undef, nup)

    @inbounds for i in 1:nup
        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        indir_up_left[i] = Vector{Int}(undef, length(sparse_row_support))
        indir_up_right[i] = Vector{Int}(undef, length(sparse_row_support))
        l, r = 1, 1
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                indir_up_left[i][l] = j
                l += 1
            else
                indir_up_right[i][r] = j
                r += 1
            end
        end
        resize!(indir_up_left[i], l - 1)
        resize!(indir_up_right[i], r - 1)
    end
    indir_low_left, indir_low_right =
        Vector{Vector{Int}}(undef, nlow), Vector{Vector{Int}}(undef, nlow)
    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        indir_low_left[i] = Vector{Int}(undef, length(sparse_row_support))
        indir_low_right[i] = Vector{Int}(undef, length(sparse_row_support))
        l, r = 1, 1
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                indir_low_left[i][l] = j
                l += 1
            else
                indir_low_right[i][r] = j
                r += 1
            end
        end
        resize!(indir_low_left[i], l - 1)
        resize!(indir_low_right[i], r - 1)
    end

    # Compute hashes of rows
    hash_vector = matrix.hash_vector
    resize!(hash_vector, nright)
    @inbounds for i in 1:nright
        hash_vector[i] = mod_p(rand(rng, AccumType), arithmetic) % CoeffType
    end
    hashes_up = Vector{AccumType}(undef, nup)
    hashes_low = Vector{AccumType}(undef, nlow)

    @inbounds for i in 1:nup
        row_hash = zero(AccumType)
        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        for j in 1:length(indir_up_right[i])
            idx = indir_up_right[i][j]
            col = sparse_row_support[idx]
            row_hash = mod_p(
                row_hash + AccumType(sparse_row_coeffs[idx]) * hash_vector[col - nleft],
                arithmetic
            )
        end
        hashes_up[i] = row_hash
    end
    @inbounds for i in 1:nlow
        row_hash = zero(AccumType)
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        for j in 1:length(indir_low_right[i])
            idx = indir_low_right[i][j]
            col = sparse_row_support[idx]
            row_hash = mod_p(
                row_hash + AccumType(sparse_row_coeffs[idx]) * hash_vector[col - nleft],
                arithmetic
            )
        end
        hashes_low[i] = row_hash
    end

    # @info "hash vector" hash_vector
    # @info "Hashes are" hashes_low hashes_up

    # Allocate the buffers
    # row = zeros(AccumType, ncols)
    row_left = zeros(AccumType, nleft)
    row_right = zeros(AccumType, ncols)
    mults = Vector{Tuple{Int, CoeffType}}(undef, nleft)

    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        hash_low = hashes_low[i]
        tmp = indir_low_left[i]
        for j in 1:length(tmp)
            idx = tmp[j]
            col = sparse_row_support[idx]
            row_left[col] = sparse_row_coeffs[idx]
        end
        # println("row_left: ", row_left)

        # for j in 1:nleft
        #     mults[j] = zero(CoeffType)
        # end
        nmults = 0
        first_nnz_column = sparse_row_support[1]
        for j in first_nnz_column:nleft
            if iszero(row_left[j])
                continue
            end

            if !isassigned(pivots, j)
                continue
            end

            # @info "reducing $i with $j"

            indices = pivots[j]
            coeffs = basis.coeffs[matrix.upper_to_coeffs[j]]

            mult = divisor(arithmetic) - row_left[j]
            nmults += 1
            mults[nmults] = (j, mult)
            hash_up = hashes_up[j]
            hash_low = mod_p(hash_low + AccumType(mult) * AccumType(hash_up), arithmetic)

            tmp = indir_up_left[j]
            for w in 1:length(tmp)
                idx = tmp[w]
                col = indices[idx]
                a = row_left[col] + AccumType(mult) * AccumType(coeffs[idx])
                row_left[col] = mod_p(a, arithmetic)
            end
        end
        # @info "reduced hash" hash_low
        # println("mults: ", mults)

        # if fully reduced
        if iszero(hash_low)
            __ZERO_HASH[] += 1
            continue
        end

        # reduce rhs by lhs pivots
        for j in 1:nright
            row_right[nleft + j] = zero(AccumType)
        end
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        tmp = indir_low_right[i]
        for j in 1:length(tmp)
            idx = tmp[j]
            col = sparse_row_support[idx]
            row_right[col] = sparse_row_coeffs[idx]
        end
        for _idx in 1:nmults
            j, mult = mults[_idx]

            indices = pivots[j]
            coeffs = basis.coeffs[matrix.upper_to_coeffs[j]]
            @invariant length(indices) == length(coeffs)

            tmp = indir_up_right[j]
            for w in 1:length(tmp)
                idx = tmp[w]
                col = indices[idx]
                row_right[col] = mod_p(
                    row_right[col] + AccumType(mult) * AccumType(coeffs[idx]),
                    arithmetic
                )
            end
        end

        # println("row_right reduced to: ", row_right)

        # Must be nonzero.

        # reduce rhs by rhs pivots
        pivot = 0
        n_nonzeros = 0
        for j in (nleft + 1):ncols
            if iszero(row_right[j])
                continue
            end

            if !isassigned(pivots, j)
                n_nonzeros += 1
                if iszero(pivot)
                    pivot = j
                end
                continue
            end

            indices = pivots[j]
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[j]]
            @invariant length(indices) == length(coeffs)

            mult = divisor(arithmetic) - row_right[j]
            for w in 1:length(indices)
                idx = indices[w]
                a = row_right[idx] + AccumType(mult) * AccumType(coeffs[w])
                row_right[idx] = mod_p(a, arithmetic)
            end
        end

        # println("after all reduction: ", row_right)

        # if fully reduced
        if n_nonzeros == 0
            __ZERO_REDUCED[] += 1
            continue
        end

        resize!(new_sparse_row_support, n_nonzeros)
        resize!(new_sparse_row_coeffs, n_nonzeros)
        linalg_extract_sparse_row!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row_right,
            convert(Int, pivot),
            ncols
        )

        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        # println(new_sparse_row_support)
        # println(new_sparse_row_coeffs)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    println("zero: ", __ZERO_REDUCED[], ", ", __ZERO_HASH[])
    __ZERO_REDUCED[] = 0
    __ZERO_HASH[] = 0

    true
end

function linalg_randomized_hashcolumns_reduce_matrix_lower_part_xxx!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nup, nlow = matrix_nrows_filled(matrix)
    nleft, nright = matrix_ncols_filled(matrix)

    # Prepare the matrix
    pivots, row_index_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # Compute hashes of rows
    hash_vector = matrix.hash_vector
    resize!(hash_vector, nright)
    @inbounds for i in 1:nright
        hash_vector[i] = mod_p(rand(rng, AccumType), arithmetic) % CoeffType
    end
    hashes_up = Vector{AccumType}(undef, nup)
    hashes_low = Vector{AccumType}(undef, nlow)
    @inbounds for i in 1:nup
        row_hash = zero(AccumType)
        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash + AccumType(sparse_row_coeffs[j]) * hash_vector[idx - nleft],
                arithmetic
            )
        end
        hashes_up[i] = row_hash
    end
    @inbounds for i in 1:nlow
        row_hash = zero(AccumType)
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash + AccumType(sparse_row_coeffs[j]) * hash_vector[idx - nleft],
                arithmetic
            )
        end
        hashes_low[i] = row_hash
    end

    # @info "hash vector" hash_vector
    # @info "Hashes are" hashes_low hashes_up

    # Allocate the buffers
    # row = zeros(AccumType, ncols)
    row_left = zeros(AccumType, nleft)
    row_right = zeros(AccumType, ncols)
    mults = Vector{Tuple{Int, CoeffType}}(undef, nleft)

    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        @invariant length(sparse_row_support) == length(sparse_row_coeffs)

        hash_low = hashes_low[i]
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx > nleft
                continue
            end
            row_left[idx] = sparse_row_coeffs[j]
        end
        # println("row_left: ", row_left)

        # for j in 1:nleft
        #     mults[j] = zero(CoeffType)
        # end
        nmults = 0
        first_nnz_column = sparse_row_support[1]
        for j in first_nnz_column:nleft
            if iszero(row_left[j])
                continue
            end

            if !isassigned(pivots, j)
                continue
            end

            # @info "reducing $i with $j"

            indices = pivots[j]
            coeffs = basis.coeffs[matrix.upper_to_coeffs[j]]

            mult = divisor(arithmetic) - row_left[j]
            nmults += 1
            mults[nmults] = (j, mult)
            hash_up = hashes_up[j]
            hash_low = mod_p(hash_low + AccumType(mult) * AccumType(hash_up), arithmetic)

            @inbounds for w in 1:length(indices)
                idx = indices[w]
                if idx > nleft
                    continue
                end
                a = row_left[idx] + AccumType(mult) * AccumType(coeffs[w])
                row_left[idx] = mod_p(a, arithmetic)
            end
        end
        # @info "reduced hash" hash_low
        # println("mults: ", mults)

        # if fully reduced
        if iszero(hash_low)
            continue
        end

        # reduce rhs by lhs pivots
        for j in 1:nright
            row_right[nleft + j] = zero(AccumType)
        end
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                continue
            end
            row_right[idx] = sparse_row_coeffs[j]
        end
        # println("row_right is: ", row_right)

        for _idx in 1:nmults
            j, mult = mults[_idx]

            indices = pivots[j]
            coeffs = basis.coeffs[matrix.upper_to_coeffs[j]]
            @invariant length(indices) == length(coeffs)

            for w in 1:length(indices)
                idx = indices[w]
                if idx <= nleft
                    continue
                end
                row_right[idx] = mod_p(
                    row_right[idx] + AccumType(mult) * AccumType(coeffs[w]),
                    arithmetic
                )
            end
        end

        # println("row_right reduced to: ", row_right)

        # Must be nonzero.

        # reduce rhs by rhs pivots
        pivot = 0
        n_nonzeros = 0
        for j in (nleft + 1):ncols
            if iszero(row_right[j])
                continue
            end

            if !isassigned(pivots, j)
                n_nonzeros += 1
                if iszero(pivot)
                    pivot = j
                end
                continue
            end

            indices = pivots[j]
            coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[j]]
            @invariant length(indices) == length(coeffs)

            mult = divisor(arithmetic) - row_right[j]
            for w in 1:length(indices)
                idx = indices[w]
                a = row_right[idx] + AccumType(mult) * AccumType(coeffs[w])
                row_right[idx] = mod_p(a, arithmetic)
            end
        end

        # println("after all reduction: ", row_right)

        # if fully reduced
        if n_nonzeros == 0
            continue
        end

        resize!(new_sparse_row_support, n_nonzeros)
        resize!(new_sparse_row_coeffs, n_nonzeros)
        linalg_extract_sparse_row!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row_right,
            convert(Int, pivot),
            ncols
        )

        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        # println(new_sparse_row_support)
        # println(new_sparse_row_coeffs)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    true
end

###

function linalg_reduce_dense_row_by_pivots_sparsedense!(
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

        linalg_vector_addmul_densedense!(row2, reducer, mult, arithmetic)
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
