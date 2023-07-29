# Macaulay matrix

# MacaulayMatrix is a sparse matrix with columns labeled by distinct monomials.
# MacaulayMatrix represents a block matrix of the following structure:
#
#   | A  B |
#   | C  D | 
# 
# The block A is in row echelon form.
# The primary operation is linear reduction of the lower part of the matrix (the
# block CD) with respect to the upper part (the block AB).

const ColumnLabel = Int32

mutable struct MacaulayMatrix{T <: Coeff}
    # A row in the matrix is represented sparsely with a vector of indices of
    # nonzero columns.
    # The coefficients of a matrix row can be stored in two different places:
    # - Either as a pointer to a polynomial in the basis that has the same
    #   coefficients.
    # - Or stored explicitly in the matrix.

    # Rows from the block AB.
    upper_rows::Vector{Vector{ColumnLabel}}
    # Rows from the block CD
    lower_rows::Vector{Vector{ColumnLabel}}

    # Explicitly stored coefficients of *some* rows
    coeffs::Vector{Vector{T}}

    # Maps column index in the range {1 ... ncols} to a corresponding monomial
    column_to_monom::Vector{MonomIdx}

    # The total number of allocated rows
    size::Int
    # The number of filled rows
    nrows::Int
    # The number of columns
    ncolumns::Int
    # The number of pivots, i.e. the number of linearly independent elements in
    # the block CD after the reduction by the block AB
    npivots::Int

    # The number of upper rows (in the AB block)
    nupper::Int
    # The number of lower rows (in the CD block)
    nlower::Int
    # The number of left columns (in the A and C blocks)
    nleft::Int
    # The number of right columns (in the B and D blocks)
    nright::Int

    # Index of the row --> position of the array of coefficients for this row
    upper_to_coeffs::Vector{Int}
    lower_to_coeffs::Vector{Int}

    # Index of the row --> monomial multiplier
    upper_to_mult::Vector{MonomIdx}
    lower_to_mult::Vector{MonomIdx}
end

# Initializes an empty matrix with coefficients of type T
function initialize_matrix(ring::PolyRing, ::Type{T}) where {T <: Coeff}
    upper_rows = Vector{Vector{ColumnLabel}}(undef, 0)
    lower_rows = Vector{Vector{ColumnLabel}}(undef, 0)
    column_to_monom = Vector{MonomIdx}(undef, 0)
    coeffs = Vector{Vector{T}}(undef, 0)
    size = 0
    npivots = 0
    nrows = 0
    ncolumns = 0
    nupper = 0
    nlower = 0
    nleft = 0
    nright = 0
    upper_to_coeffs = Vector{Int}(undef, 0)
    lower_to_coeffs = Vector{Int}(undef, 0)
    upper_to_mult = Vector{MonomIdx}(undef, 0)
    lower_to_mult = Vector{MonomIdx}(undef, 0)
    MacaulayMatrix(
        upper_rows,
        lower_rows,
        coeffs,
        column_to_monom,
        size,
        nrows,
        ncolumns,
        npivots,
        nupper,
        nlower,
        nleft,
        nright,
        upper_to_coeffs,
        lower_to_coeffs,
        upper_to_mult,
        lower_to_mult
    )
end

# Returns a string with human-readable information about the matrix. Can be used
# in logging
function repr_matrix(matrix::MacaulayMatrix{T}) where {T}
    m, n = matrix.nrows, matrix.ncolumns
    m_A, n_A = matrix.nupper, matrix.nleft
    m_B, n_B = matrix.nupper, matrix.nright
    m_C, n_C = matrix.nlower, matrix.nleft
    m_D, n_D = matrix.nlower, matrix.nright
    nnz_A, nnz_B, nnz_C, nnz_D = 0, 0, 0, 0
    A_ref, A_rref = true, true
    max_canvas_width = 40
    canvas = CanvasMatrix2x2(m_A, m_C, n_A, n_B, max_width=max_canvas_width)
    for i in 1:(matrix.nupper)
        row = matrix.upper_rows[i]
        if length(row) > 1
            A_rref = false
        end
        if row[1] < i
            A_ref = true
        end
        for j in 1:length(row)
            point!(canvas, i, row[j])
            if row[j] <= matrix.nleft
                nnz_A += 1
            else
                nnz_B += 1
            end
        end
    end
    for i in 1:(matrix.nlower)
        row = matrix.lower_rows[i]
        for j in 1:length(row)
            point!(canvas, matrix.nupper + i, row[j])
            if row[j] <= matrix.nleft
                nnz_C += 1
            else
                nnz_D += 1
            end
        end
    end
    nnz = nnz_A + nnz_B + nnz_C + nnz_D
    percent(x) = round(100 * x, digits=2)
    s = """
    $(typeof(matrix))
    $m x $n with $nnz nonzeros ($(percent(nnz / (m * n))) %)
    A: $(m_A) x $(n_A) with $(nnz_A) nonzeros (REF: $(A_ref), RREF: $(A_rref))
    B: $(m_B) x $(n_B) with $(nnz_B) nonzeros
    C: $(m_C) x $(n_C) with $(nnz_C) nonzeros
    D: $(m_D) x $(n_D) with $(nnz_D) nonzeros

    Sparsity pattern:

    $(canvas)
    """
    s
end

# Checks that the matrix is well formed.
function matrix_well_formed(key, matrix::MacaulayMatrix)
    if !(matrix.size == length(matrix.upper_rows))
        return false
    end
    if !(matrix.ncolumns == matrix.nleft + matrix.nright)
        return false
    end
    if !(matrix.nrows == matrix.nupperper + matrix.nlower)
        return false
    end
    if key in (:in_reduction_apply!, :in_reduction!, :in_reduction_learn!)
        for i in 1:length(matrix.upper_rows)
            !isassigned(matrix.upper_rows, i) && continue
            if isempty(matrix.upper_rows[i])
                return false
            end
        end
        for i in 1:length(matrix.lower_rows)
            !isassigned(matrix.lower_rows, i) && continue
            if isempty(matrix.lower_rows[i])
                return false
            end
        end
    end
    true
end

# Refresh and initialize matrix
function reinitialize_matrix!(matrix::MacaulayMatrix{T}, size::Int) where {T}
    # TODO: do not reinitialize when the present size is enough 
    resize!(matrix.upper_rows, size * 2)
    resize!(matrix.lower_rows, size * 2)
    resize!(matrix.upper_to_coeffs, size * 2)
    resize!(matrix.lower_to_coeffs, size * 2)
    resize!(matrix.upper_to_mult, size * 2)
    resize!(matrix.lower_to_mult, size * 2)
    matrix.size = 2 * size
    matrix.ncolumns = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nupper = 0
    matrix.nlower = 0
    matrix
end

# Do linear reduction of the lower part of the matrix
function linear_algebra!(
    ring::PolyRing,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::Symbol,
    rng
)
    resize!(matrix.coeffs, matrix.nlower)
    # @invariant matrix_well_formed(:none, matrix)
    if linalg === :deterministic
        deterministic_sparse_rref!(ring, matrix, basis)
        return true
    else
        @assert linalg === :randomized
        randomized_sparse_rref!(ring, matrix, basis, rng)
        return true
    end
    true
end

# Do linear reduction of the lower part of the matrix.
# Returns true if the reduction is successful
function linear_algebra!(
    graph::ComputationGraphF4,
    ring::PolyRing,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::Symbol,
    rng
)
    resize!(matrix.coeffs, matrix.nlower)
    # @invariant matrix_well_formed(:none, matrix)
    if linalg === :learn
        return learn_sparse_rref!(graph, ring, matrix, basis, rng)
    else
        @assert linalg === :apply
        return apply_sparse_rref!(graph, ring, matrix, basis, rng)
    end
    true
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_sparse_row!(
    row::Vector{T},
    arithmetic::A
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = mod_x(invmod(row[1], divisor(arithmetic)), arithmetic)
    @inbounds row[1] = one(row[1])
    @inbounds for i in 2:length(row)
        row[i] = mod_x(row[i] * pinv, arithmetic)
    end
    row
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_sparse_row!(
    row::Vector{T},
    arithmetic::A
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = inv(row[1])
    @inbounds row[1] = one(row[1])
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end
    row
end

function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{ColumnLabel},
    coeffs::Vector{T},
    arithmetic::A
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @invariant isone(coeffs[1])
    @inbounds mul = divisor(arithmetic) - row[indices[1]]

    # On our benchmarks, usually,
    #   length(row) / length(indices)
    # roughly varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = mod_x(row[idx] + mul * coeffs[j], arithmetic)
    end

    nothing
end

function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{ColumnLabel},
    coeffs::Vector{T},
    arithmetic::A
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @invariant isone(coeffs[1])

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
        # Base.GMP.MPQ.mul!(buf2, cfs[j], mul)
        # row[idx] = row[idx] + buf2
    end

    nothing
end

# Load coefficients from `coeffs` into dense `row` at `indices`
function load_sparse_row!(row::Vector{T}, indices, coeffs) where {T <: CoeffFF}
    @inbounds for i in 1:length(row)
        row[i] = T(0)
    end
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end
end

# Load coefficients from `coeffs` into dense `row` at `indices`
function load_sparse_row!(row::Vector{T}, indices, coeffs) where {T <: CoeffQQ}
    row .= T(0)
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end
end

# Traverses the dense `row` at positions `from..to` and extracts all nonzero
# entries to the given sparse row  
function extract_sparse_row!(indices, coeffs, row::Vector{T}, from::Int, to::Int) where {T}
    # NOTE: assumes the sparse row has the necessary capacity
    z = zero(T)
    j = 1
    @inbounds for i in from:to
        if row[i] != z
            indices[j] = i
            coeffs[j] = row[i]
            j += 1
        end
    end
    nothing
end

# function extract_sparse_row_2!(vals::Vector{UInt}, dense::Vector{UInt})
#     j = 1
#     z = Vec{4, UInt}((0, 0, 0, 0))
#     @inbounds for i in 1:8:length(dense)
#         v1 = SIMD.vload(Vec{4, UInt}, dense, i)
#         v2 = SIMD.vload(Vec{4, UInt}, dense, i + 4)
#         mask1 = !iszero(v1)
#         mask2 = !iszero(v2)
#         j1 = sum(mask1)
#         j2 = sum(mask2)
#         vstorec(v1, vals, j, mask1)
#         vstorec(v2, vals, j + j1, mask2)
#         j += j1 + j2
#     end
#     nothing
# end

function reduce_dense_row_by_pivots_sparse!(
    row::Vector{T},
    matrix::MacaulayMatrix{T},
    basis::Basis{T},
    pivots::Vector{Vector{ColumnLabel}},
    start_column::ColumnLabel,
    tmp_pos::ColumnLabel,
    arithmetic::A;
    exact_column_mapping=false
) where {T <: Coeff, A <: AbstractArithmetic}
    ncolumns = matrix.ncolumns
    nleft = matrix.nleft
    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero elemen in the reduced row
    new_pivot = -1

    # TODO
    @inbounds for i in start_column:ncolumns
        # @inbounds for i in start_column:nleft
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
        else
            # if reducer is from the lower part of the matrix
            @log level = -777 "Lower part" i
            coeffs = matrix.coeffs[matrix.lower_to_coeffs[i]]
        end

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)
        @invariant iszero(row[i])
    end
    # form and return the resulting row in sparse format 
    indices = Vector{ColumnLabel}(undef, nonzeros)
    coeffs = Vector{T}(undef, nonzeros)
    # all reduced to zero!
    if nonzeros == 0
        return true, indices, coeffs
    end
    extract_sparse_row!(indices, coeffs, row, Int(start_column), Int(ncolumns))
    return false, indices, coeffs
end

# same as above, but remembers the rows that acted as reducers and stores them
# in active_reducers
function reduce_dense_row_by_pivots_sparse!(
    active_reducers,
    row::Vector{T},
    matrix::MacaulayMatrix{T},
    basis::Basis{T},
    pivots::Vector{Vector{ColumnLabel}},
    start_column::ColumnLabel,
    tmp_pos::ColumnLabel,
    arithmetic::A;
    exact_column_mapping=false
) where {T <: Coeff, A <: AbstractArithmetic}
    ncolumns = matrix.ncolumns
    nleft = matrix.nleft
    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero elemen in the reduced row
    new_pivot = -1

    @inbounds for i in start_column:ncolumns
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            push!(active_reducers, (matrix.upper_to_coeffs[i], matrix.upper_to_mult[i]))
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.coeffs[matrix.lower_to_coeffs[i]]
        end

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)
        @invariant iszero(row[i])
    end
    # form and return the resulting row in sparse format 
    indices = Vector{ColumnLabel}(undef, nonzeros)
    coeffs = Vector{T}(undef, nonzeros)
    # all reduced to zero!
    if nonzeros == 0
        return true, indices, coeffs
    end
    extract_sparse_row!(indices, coeffs, row, Int(start_column), Int(ncolumns))
    return false, indices, coeffs
end

function get_absolute_pivots!(matrix::MacaulayMatrix)
    pivots = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nupper)
        pivots[i] = matrix.upper_rows[i]
        @invariant pivots[i][1] == i
    end
    l2c_tmp = Vector{ColumnLabel}(undef, max(matrix.ncolumns, matrix.nlower))
    @inbounds for i in 1:(matrix.nlower)
        l2c_tmp[matrix.lower_rows[i][1]] = matrix.lower_to_coeffs[i]
    end
    pivots, l2c_tmp
end

function interreduce_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{ColumnLabel}},
    arithmetic;
    reversed_rows=false
) where {C <: Coeff}
    # the number of new pivots
    nnew_pivots = 0
    row = zeros(C, matrix.ncolumns)
    resize!(matrix.lower_rows, matrix.nright)
    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:(matrix.nright)
        k = matrix.ncolumns - i + 1
        !isassigned(pivots, k) && continue

        indices = pivots[k]
        if k <= matrix.nleft
            # upper part of matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[k]]
        else
            # lower part of matrix
            coeffs = matrix.coeffs[matrix.lower_to_coeffs[k]]
        end
        @assert length(coeffs) == length(indices)

        start_column = indices[1]
        load_sparse_row!(row, indices, coeffs)

        _, new_indices, new_coeffs = reduce_dense_row_by_pivots_sparse!(
            row,
            matrix,
            basis,
            pivots,
            ColumnLabel(start_column),
            ColumnLabel(start_column),
            arithmetic
        )
        nnew_pivots += 1

        # update row and coeffs
        if !reversed_rows
            matrix.lower_rows[nnew_pivots] = new_indices
            matrix.coeffs[matrix.lower_to_coeffs[k]] = new_coeffs
            pivots[k] = matrix.lower_rows[nnew_pivots]
        else
            matrix.lower_rows[matrix.nrows - nnew_pivots + 1] = new_indices
            matrix.coeffs[matrix.lower_to_coeffs[k]] = new_coeffs
            pivots[k] = matrix.lower_rows[matrix.nrows - nnew_pivots + 1]
        end
    end

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = nnew_pivots
    resize!(matrix.lower_rows, nnew_pivots)
end

function interreduce_lower_part_learn!(
    context::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{ColumnLabel}},
    arithmetic;
    reversed_rows=false
) where {C <: Coeff}
    # the number of new pivots
    nnew_pivots = 0
    row = zeros(C, matrix.ncolumns)
    resize!(matrix.lower_rows, matrix.nright)
    # Indices of rows that did no reduce to zero
    not_reduced_to_zero = Vector{Int}()
    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:(matrix.nright)
        k = matrix.ncolumns - i + 1
        !isassigned(pivots, k) && continue

        indices = pivots[k]
        if k <= matrix.nleft
            # upper part of matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[k]]
        else
            # lower part of matrix
            coeffs = matrix.coeffs[matrix.lower_to_coeffs[k]]
        end
        @assert length(coeffs) == length(indices)

        start_column = indices[1]
        load_sparse_row!(row, indices, coeffs)

        _, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            row,
            matrix,
            basis,
            pivots,
            ColumnLabel(start_column),
            ColumnLabel(start_column),
            arithmetic
        )
        nnew_pivots += 1

        push!(not_reduced_to_zero, k)

        # update row and coeffs
        if !reversed_rows
            matrix.lower_rows[nnew_pivots] = newrow
            matrix.coeffs[matrix.lower_to_coeffs[k]] = newcfs
            pivots[k] = matrix.lower_rows[nnew_pivots]
        else
            matrix.lower_rows[matrix.nrows - nnew_pivots + 1] = newrow
            matrix.coeffs[matrix.lower_to_coeffs[k]] = newcfs
            pivots[k] = matrix.lower_rows[matrix.nrows - nnew_pivots + 1]
        end
    end

    # Update the context
    push!(
        context.matrix_infos,
        (nup=matrix.nupper, nlow=matrix.nlower, ncols=matrix.ncolumns)
    )
    push!(context.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        context.matrix_upper_rows,
        (
            matrix.upper_to_coeffs[1:(matrix.nupper)],
            matrix.upper_to_mult[1:(matrix.nupper)]
        )
    )
    push!(context.matrix_lower_rows, (Int[], Int[]))

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = nnew_pivots
    resize!(matrix.lower_rows, nnew_pivots)
end

# Linear algebra option 1
function deterministic_sparse_rref!(
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lower_rows
    densecoeffs = zeros(C, ncols)

    @log level = -777 "" matrix.nupper matrix.nlower matrix.nleft matrix.nright
    @log level = -777 "" matrix.ncolumns matrix.nrows matrix.size
    @log level = -777 "" map(first, pivs[1:(matrix.nupper)]) == collect(1:(matrix.nupper))

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]

        # we load coefficients into dense array
        # into rowexps indices
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # if fully reduced
        zeroed && continue

        @log level = -777 "" newrow[1]

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
end

function learn_sparse_rref!(
    graph::ComputationGraphF4,
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    rng
) where {C <: Coeff}
    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    @log level = -6 "absolute_pivots!" pivs l2c_tmp

    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lower_rows
    densecoeffs = zeros(C, ncols)

    not_reduced_to_zero = Vector{Int}()
    useful_reducers = Set{Tuple{Int, MonomIdx}}()

    @log level = -6 "Low to coef" rowidx2coef matrix.lower_to_mult

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]

        @log level = -7 "$i from $nlow" rowexps cfsref rowidx2coef

        # we load coefficients into dense array
        # into rowexps indices
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        reducers = Any[]
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            reducers,
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # if fully reduced
        zeroed && continue

        # NOTE: we are not adding reducers from lowrows!

        push!(not_reduced_to_zero, i)
        for rr in reducers
            push!(useful_reducers, rr)
        end

        @log level = -7 "Not zero" i rowidx2coef[i] matrix.lower_to_mult[i]

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    useful_reducers_sorted = sort(collect(useful_reducers))
    @log level = -7 "" useful_reducers_sorted
    push!(
        graph.matrix_infos,
        (nup=matrix.nupper, nlow=matrix.nlower, ncols=matrix.ncolumns)
    )
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (map(first, useful_reducers_sorted), map(last, useful_reducers_sorted))
    )
    push!(
        graph.matrix_lower_rows,
        (rowidx2coef[not_reduced_to_zero], matrix.lower_to_mult[not_reduced_to_zero])
    )

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
    true
end

function apply_sparse_rref!(
    graph::ComputationGraphF4,
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    rng
) where {C <: Coeff}
    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    @log level = -7 "absolute_pivots!" pivs l2c_tmp

    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lower_rows
    densecoeffs = zeros(C, ncols)

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]

        # we load coefficients into dense array
        # into rowexps indices
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # if fully reduced
        if zeroed
            # then something has gone wrong in tracing 
            return false
        end

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
    true
end

function nblocks_in_randomized(nrows::Int)
    floor(Int, sqrt(nrows / 3)) + 1
end

function randomized_sparse_rref!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    rng
) where {C <: Coeff}
    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lower_rows

    nblocks = nblocks_in_randomized(nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem

    densecoeffs = zeros(C, ncols)
    mulcoeffs = zeros(C, rowsperblock)

    for i in 1:nblocks
        nrowsupper = nlow > i * rowsperblock ? i * rowsperblock : nlow
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue # end right here

        ctr = 0
        @inbounds while ctr < nrowstotal
            for j in 1:nrowstotal
                mulcoeffs[j] = mod_x(rand(rng, C), arithmetic)
            end
            densecoeffs .= C(0)
            startcol = ncols

            @inbounds for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                rowexps = upivs[rowidx]
                cfsref = basis.coeffs[rowidx2coef[rowidx]]
                startcol = min(startcol, rowexps[1])

                for l in 1:length(rowexps)
                    ridx = rowexps[l]
                    densecoeffs[ridx] =
                        mod_x(densecoeffs[ridx] + mulcoeffs[k] * cfsref[l], arithmetic)
                end
            end

            # reduce it with known pivots from matrix.upper_rows
            # first nonzero in densecoeffs is at startcol position
            zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
                densecoeffs,
                matrix,
                basis,
                pivs,
                ColumnLabel(startcol),
                ColumnLabel(-1),
                arithmetic
            )

            if zeroed
                ctr = nrowstotal
                break
            end

            absolute_i = (i - 1) * rowsperblock + ctr + 1

            # matrix coeffs sparsely stores coefficients of new row
            matrix.coeffs[absolute_i] = newcfs
            # add new pivot at column index newrow[1]
            #  (which is the first nnz column of newrow)
            pivs[newrow[1]] = newrow
            # set ref to coefficient to matrix
            # guaranteed to be from lower part
            matrix.lower_to_coeffs[newrow[1]] = absolute_i

            # normalize if needed
            normalize_sparse_row!(matrix.coeffs[absolute_i], arithmetic)

            ctr += 1
        end
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
end

function deterministic_sparse_rref_interreduce_apply!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lower_rows, matrix.ncolumns)
    resize!(matrix.upper_to_coeffs, matrix.ncolumns)
    resize!(matrix.upper_to_mult, matrix.ncolumns)
    resize!(matrix.lower_to_coeffs, matrix.ncolumns)
    resize!(matrix.lower_to_mult, matrix.ncolumns)
    resize!(matrix.coeffs, matrix.ncolumns)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
        matrix.lower_to_coeffs[matrix.upper_rows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.upper_to_coeffs[i]])
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic, reversed_rows=true)
end

function deterministic_sparse_rref_interreduce!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lower_rows, matrix.ncolumns)
    resize!(matrix.upper_to_coeffs, matrix.ncolumns)
    resize!(matrix.upper_to_mult, matrix.ncolumns)
    resize!(matrix.lower_to_coeffs, matrix.ncolumns)
    resize!(matrix.lower_to_mult, matrix.ncolumns)
    resize!(matrix.coeffs, matrix.ncolumns)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
        matrix.lower_to_coeffs[matrix.upper_rows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.upper_to_coeffs[i]])
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic, reversed_rows=true)
    true
end

function deterministic_sparse_rref_interreduce_learn!(
    graph,
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lower_rows, matrix.ncolumns)
    resize!(matrix.upper_to_coeffs, matrix.ncolumns)
    resize!(matrix.upper_to_mult, matrix.ncolumns)
    resize!(matrix.lower_to_coeffs, matrix.ncolumns)
    resize!(matrix.lower_to_mult, matrix.ncolumns)
    resize!(matrix.coeffs, matrix.ncolumns)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
        matrix.lower_to_coeffs[matrix.upper_rows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.upper_to_coeffs[i]])
    end

    interreduce_lower_part_learn!(
        graph,
        matrix,
        basis,
        pivs,
        arithmetic,
        reversed_rows=true
    )
end

function deterministic_sparse_rref_isgroebner!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    upivs = matrix.lower_rows

    densecoeffs = zeros(C, matrix.ncolumns)

    @inbounds for i in 1:(matrix.nlower)
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]
        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, _, _ = reduce_dense_row_by_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # if fully reduced
        zeroed && continue
        return false
    end
    return true
end

function deterministic_sparse_rref_nf!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    tobereduced::Basis{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.coeffs, matrix.nlower)

    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    upivs = matrix.lower_rows

    densecoeffs = zeros(C, ncols)

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = tobereduced.coeffs[rowidx2coef[i]]
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # # @warn "reduced " zeroed
        # if fully reduced
        zeroed && continue

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        matrix.lower_rows[i] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[i] = i
    end
    matrix.npivots = matrix.nrows = matrix.size = matrix.nlower
end

function column_to_monom_mapping!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    column_to_monom = Vector{MonomIdx}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        # column to hash index
        column_to_monom[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == 2
            k += 1
        end
    end

    sort_columns_by_labels!(column_to_monom, symbol_ht)

    matrix.nleft = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        hdata[column_to_monom[k]].idx = k
    end

    @inbounds for k in 1:(matrix.nupper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nlower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.ncolumns = matrix.nleft + matrix.nright

    @assert matrix.nleft + matrix.nright == symbol_ht.load - 1 == matrix.ncolumns
    @assert matrix.nlower + matrix.nupper == matrix.nrows

    matrix.column_to_monom = column_to_monom
end

function convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)

    # we mutate basis array directly by adding new elements

    resize_basis_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.nprocessed

    @inbounds for i in 1:(matrix.npivots)
        colidx = rows[i][1]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.column_to_monom)
        basis.coeffs[crs + i] = matrix.coeffs[matrix.lower_to_coeffs[colidx]]
        basis.monoms[crs + i] = matrix.lower_rows[i]
    end

    basis.nfilled += matrix.npivots
end

function convert_rows_to_basis_elements_nf!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    resize_basis_if_needed!(basis, matrix.npivots)

    @inbounds for i in 1:(matrix.npivots)
        basis.nprocessed += 1
        basis.nnonredundant += 1
        basis.nonredundant[basis.nnonredundant] = basis.nprocessed
        if isassigned(matrix.coeffs, i)
            row = matrix.lower_rows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[basis.nprocessed] = matrix.coeffs[i]
            basis.monoms[basis.nprocessed] = row
        else
            empty!(basis.coeffs[basis.nprocessed])
            empty!(basis.monoms[basis.nprocessed])
        end
    end

    nothing
end

function insert_in_basis_hash_table_pivots(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomIdx}
) where {M}
    resize_hashtable_if_needed!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.monoms
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = column_to_monom[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = next_lookup_index(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if ishashcollision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask, sdata[hidx].deg)

        ht.load += 1
    end

    nothing
end
