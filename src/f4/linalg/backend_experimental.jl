# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

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
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_direct_rref_sparse!"
    @log :matrix matrix_string_repr(matrix)

    # Produce the RREF of AB
    linalg_interreduce_matrix_upper_part!(matrix, basis, arithmetic)
    # Use the produced RREF to perform reduction of CD
    linalg_reduce_matrix_lower_part!(matrix, basis, arithmetic)
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
    resize!(matrix.B_coeffs_dense, ncols)

    # Compute hashes of rows in the B block
    hash_vector = matrix.buffer_hash_vector
    resize!(hash_vector, nright)
    @inbounds for i in 1:nright
        hash_vector[i] = mod_p(rand(AccumType), arithmetic)
    end
    @inbounds for i in 1:nup
        row_hash = zero(AccumType)
        sparse_row_support = matrix.upper_rows[i]
        sparse_row_coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash +
                AccumType(sparse_row_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                arithmetic
            )
        end
        matrix.B_coeffs_dense[sparse_row_support[1]] = [row_hash % CoeffType]
    end

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    row1 = zeros(AccumType, nleft)
    row2 = zeros(AccumType, 1)

    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        for j in 1:nleft
            row1[j] = zero(AccumType)
        end
        row_hash = zero(AccumType)
        for j in 1:length(sparse_row_support)
            idx = sparse_row_support[j]
            if idx <= nleft
                row1[idx] = sparse_row_coeffs[j]
            else
                row_hash = mod_p(
                    row_hash +
                    AccumType(sparse_row_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                    arithmetic
                )
            end
        end
        row2[1] = row_hash

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        first_nnz_column = sparse_row_support[1]
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

            linalg_vector_addmul_sparsedense_mod_p!(row, indices, coeffs, arithmetic)

            linalg_vector_addmul_densedense!(row2, reducer, mult, arithmetic)
        end
        zeroed = iszero(row2[1])

        # if fully reduced
        if zeroed && iszero(pivot)
            continue
        end

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # otherwise, reduce once again
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
            continue
        end

        linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        row_hash = zero(AccumType)
        for j in 1:length(new_sparse_row_support)
            idx = new_sparse_row_support[j]
            if idx <= nleft
                continue
            end
            row_hash = mod_p(
                row_hash +
                AccumType(new_sparse_row_coeffs[j]) * AccumType(hash_vector[idx - nleft]),
                arithmetic
            )
        end
        matrix.B_coeffs_dense[new_sparse_row_support[1]] = [row_hash]

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
