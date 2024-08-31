# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# High level

function linalg_deterministic_sparse_changematrix!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part

    @log :matrix "linalg_deterministic_sparse_changematrix!"
    @log :matrix matrix_string_repr(matrix)

    matrix_changematrix_initialize!(matrix, matrix.nrows_filled_lower)

    # Reduce CD with AB
    linalg_reduce_matrix_lower_part_changematrix!(matrix, basis, arithmetic)
    # Interreduce CD
    linalg_interreduce_matrix_pivots_changematrix!(matrix, basis, arithmetic)
    true
end

function linalg_deterministic_sparse_interreduction_changematrix!(
    matrix::MacaulayMatrix,
    basis::Basis{C},
    arithmetic::AbstractArithmetic
) where {C}
    sort_matrix_upper_rows!(matrix)
    @log :matrix "linalg_deterministic_sparse_interreduction_changematrix!"
    @log :matrix matrix_string_repr(matrix)

    matrix_changematrix_initialize!(matrix, matrix.nrows_filled_upper)
    for i in 1:(matrix.nrows_filled_upper)
        poly_id = matrix.upper_to_coeffs[i]
        mult_id = matrix.upper_to_mult[i]
        matrix.changematrix[i] = Dict{Tuple{Int, MonomId}, C}((poly_id, mult_id) => one(C))
    end

    # Prepare the matrix
    linalg_prepare_matrix_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    linalg_interreduce_matrix_pivots_changematrix!(
        matrix,
        basis,
        arithmetic,
        reversed_rows=true
    )
    true
end

###
# Low level

function linalg_reduce_matrix_lower_part_changematrix!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic::AbstractArithmetic{AccumType, CoeffType}
) where {CoeffType <: Coeff, AccumType <: Coeff}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)
    nup, nlow = matrix_nrows_filled(matrix)

    # Prepare the matrix
    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)
    changematrix = matrix.changematrix

    # Allocate the buffers
    reducer_rows = Vector{Tuple{Int, CoeffType}}()
    row = zeros(AccumType, ncols)
    new_sparse_row_support, new_sparse_row_coeffs = linalg_new_empty_sparse_row(CoeffType)

    @inbounds for i in 1:nlow
        sparse_row_support = matrix.lower_rows[i]
        sparse_row_coeffs = basis.coeffs[row_idx_to_coeffs[i]]

        linalg_load_sparse_row!(row, sparse_row_support, sparse_row_coeffs)

        # Additionally record the indices of rows that participated in reduction
        # of the given row
        empty!(reducer_rows)
        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse_changematrix!(
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
            
        )

        # if fully reduced
        zeroed && continue

        # NOTE: we are not recording reducers from the lower part of the matrix
        cm = changematrix[i]
        cm[(row_idx_to_coeffs[i], matrix.lower_to_mult[i])] = one(CoeffType)
        for (idx, quo) in reducer_rows
            if idx <= nleft
                poly_idx, poly_mult = matrix.upper_to_coeffs[idx], matrix.upper_to_mult[idx]
                if !haskey(cm, (poly_idx, poly_mult))
                    cm[(poly_idx, poly_mult)] = zero(CoeffType)
                end
                cm[(poly_idx, poly_mult)] =
                    mod_p(cm[(poly_idx, poly_mult)] + AccumType(quo), arithmetic)
            else
                row_idx = matrix.lower_to_coeffs[idx]
                for ((poly_idx, poly_mult), cf) in changematrix[row_idx]
                    if !haskey(cm, (poly_idx, poly_mult))
                        cm[(poly_idx, poly_mult)] = zero(CoeffType)
                    end
                    cm[(poly_idx, poly_mult)] =
                        mod_p(cm[(poly_idx, poly_mult)] + quo * AccumType(cf), arithmetic)
                end
            end
        end

        @invariant length(new_sparse_row_support) == length(new_sparse_row_coeffs)
        pinv = linalg_row_make_monic!(new_sparse_row_coeffs, arithmetic)

        cm2 = empty(cm)
        for ((poly_idx, quo_idx), quo_cf) in cm
            cm2[(poly_idx, quo_idx)] = mod_p(AccumType(pinv) * quo_cf, arithmetic)
        end
        changematrix[i] = cm2

        matrix.some_coeffs[i] = new_sparse_row_coeffs
        pivots[new_sparse_row_support[1]] = new_sparse_row_support
        matrix.lower_to_coeffs[new_sparse_row_support[1]] = i

        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
    end

    true
end

function linalg_interreduce_matrix_pivots_changematrix!(
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
    changematrix = matrix.changematrix
    compact_changematrix = similar(changematrix)

    # Allocate the buffers
    reducer_rows = Vector{Tuple{Int, CoeffType}}()
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

        empty!(reducer_rows)
        new_sparse_row_support, new_sparse_row_coeffs =
            linalg_new_empty_sparse_row(CoeffType)
        first_nnz_column = sparse_row_support[1]
        zeroed = linalg_reduce_dense_row_by_pivots_sparse_changematrix!(
            new_sparse_row_support,
            new_sparse_row_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            ncols,
            arithmetic,
            ignore_column=first_nnz_column,
            reducer_rows
        )

        # If the row reduced to zero
        if zeroed
            any_zeroed = true
            continue
        end

        current_row_idx = matrix.lower_to_coeffs[abs_column_idx]
        cm = changematrix[current_row_idx]
        for (idx, quo) in reducer_rows
            if idx <= nleft
                @assert false
                poly_idx, poly_mult = matrix.upper_to_coeffs[idx], matrix.upper_to_mult[idx]
                cm[(poly_idx, poly_mult)] = quo
            else
                row_idx = matrix.lower_to_coeffs[idx]
                for ((poly_idx, poly_mult), cf) in changematrix[row_idx]
                    if !haskey(cm, (poly_idx, poly_mult))
                        cm[(poly_idx, poly_mult)] = zero(CoeffType)
                    end
                    cm[(poly_idx, poly_mult)] =
                        mod_p(cm[(poly_idx, poly_mult)] + quo * AccumType(cf), arithmetic)
                end
            end
        end

        new_pivots += 1
        not_reduced_to_zero[new_pivots] = i

        # Update the row support and coefficients.
        # TODO: maybe get rid of the reversed_rows parameter?
        if !reversed_rows
            compact_changematrix[new_pivots] = cm
            matrix.lower_rows[new_pivots] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] =
                new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[new_pivots]
        else
            compact_changematrix[nupper - new_pivots + 1] = cm
            matrix.lower_rows[nupper - new_pivots + 1] = new_sparse_row_support
            matrix.some_coeffs[matrix.lower_to_coeffs[abs_column_idx]] =
                new_sparse_row_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[nupper - new_pivots + 1]
        end
    end

    matrix.npivots = new_pivots
    matrix.changematrix = compact_changematrix
    resize!(matrix.changematrix, new_pivots)
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)

    true, any_zeroed, not_reduced_to_zero
end

function linalg_reduce_dense_row_by_pivots_sparse_changematrix!(
    new_sparse_row_support::Vector{I},
    new_sparse_row_coeffs::Vector{C},
    row::Vector{A},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{I}},
    start_column::Integer,
    end_column::Integer,
    arithmetic::AbstractArithmeticZp{A, C},
    active_reducers;
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
        else
            pivot_coeffs = matrix.some_coeffs[matrix.lower_to_coeffs[i]]
        end
        @invariant length(pivot_support) == length(pivot_coeffs)

        mul = linalg_vector_addmul_sparsedense_mod_p!(
            row,
            pivot_support,
            pivot_coeffs,
            arithmetic
        )

        push!(active_reducers, (i, mul))

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
