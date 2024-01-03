# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# High level

function linalg_randomized_sparse!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "linalg_randomized_sparse!"
    @log level = -3 matrix_string_repr(matrix)
    # Reduce CD with AB
    linalg_randomized_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    # Interreduce CD
    linalg_interreduce_matrix_pivots!(matrix, basis, arithmetic)
    true
end

###
# Low level

# The number of blocks to split the matrix into
function linalg_nblocks_in_randomized(nrows::Int)
    floor(Int, sqrt(nrows / 3)) + 1
end

# Given a matrix of the following form, where A is in REF,
#   A B
#   C D
# reduces the lower part CD with respect to the pivots in the upper part AB.
# As a result, a matrix of the following form is produced:
#   A B
#   0 D'
function linalg_randomized_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType},
    rng::AbstractRNG
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    _, nlow = matrix_nrows_filled(matrix)
    if nlow <= 2
        @log level = -4 """
        Too few rows in the matrix. 
        Consider switching to another backend to avoid the overhead of randomization."""
    end

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, matrix.nrows_filled_lower)

    # Set up the blocks
    nblocks = linalg_nblocks_in_randomized(nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem
    @log level = -3 """
    Rows in the lower part: $nlow
    The bumber of blocks: $nblocks
    Rows per block: $rowsperblock"""

    # Allocate the buffers
    row = zeros(AccumType, ncols)
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)
    rows_multipliers = zeros(AccumType, rowsperblock)

    # Split the rows in the lower part of the matrix into several blocks. 
    # Handle each block separately. 
    # For each block, produce a random linear combination of the rows in the
    # block and reduce it by pivot rows. Repeat while nonzero rows are emerging.
    for i in 1:nblocks
        nrowsupper = nlow > i * rowsperblock ? i * rowsperblock : nlow
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        new_pivots_count = 0
        @inbounds while new_pivots_count < nrowstotal
            # Produce a random linear combination of rows in this block
            # TODO: does not work for the rationals
            # TODO: can vectorize random number generation
            # TODO: sample from a smaller range of values
            for j in 1:nrowstotal
                rows_multipliers[j] = mod_p(rand(rng, AccumType), arithmetic)
            end
            row .= AccumType(0)
            first_nnz_col = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                sparse_row_support = matrix.lower_rows[rowidx]
                sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[rowidx]]
                @invariant length(sparse_row_support) == length(sparse_row_coeffs)

                first_nnz_col = min(first_nnz_col, sparse_row_support[1])
                for l in 1:length(sparse_row_coeffs)
                    colidx = sparse_row_support[l]
                    row[colidx] = mod_p(
                        row[colidx] +
                        AccumType(rows_multipliers[k]) * AccumType(sparse_row_coeffs[l]),
                        arithmetic
                    )
                end
            end

            # Reduce the random linear combination by pivots
            @invariant 1 <= first_nnz_col <= length(row)
            zeroed = linalg_reduce_dense_row_by_pivots_sparse!(
                new_sparse_row_support,
                new_sparse_row_coeffs,
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

            linalg_normalize_row!(new_sparse_row_coeffs, arithmetic)
            @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)

            pivots[new_sparse_row_support[1]] = new_sparse_row_support

            absolute_row_index = (i - 1) * rowsperblock + new_pivots_count
            matrix.some_coeffs[absolute_row_index] = new_sparse_row_coeffs
            matrix.lower_to_coeffs[new_sparse_row_support[1]] = absolute_row_index

            new_sparse_row_support, new_sparse_row_coeffs =
                linalg_new_empty_sparse_row(CoeffType)
        end
    end

    true
end
