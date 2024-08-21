# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Macaulay matrix

# MacaulayMatrix is a sparse matrix with columns labeled by distinct monomials.
# MacaulayMatrix represents a block matrix of the following structure:
#
#   | A  B |
#   | C  D | 
#
# The sparsity structure of the matrix upon creation may look similar to:
#
#           A          B
#     ......... | ......
#       ....... | ......
#         ..... | ......
#  A        ... | ......
#             . |  .....
#     ------------------
#            .. | ......
#  C      ..... | ......
#      ........ | ......
#
# NOTE: Check out `matrix_string_repr` for printing this nice matrix representation in
# runtime.

# The primary action on a MacaulayMatrix is computing
#   D - C inv(A) B,
# which is equivalent to reducing the CD block (lower part of the matrix) with
# the AB block (upper part of the matrix) via Gaussian elimination.
#
# Usually, by construction,
# - Block A is already in row-echelon-form.
# - Block A is the largest block, and is very-very sparse.
# - Permuting columns is not allowed.

# Each column in the matrix is associated with a single integer of this type.
const ColumnLabel = Int32

# Tags for matrix columns from different blocks.
const NON_PIVOT_COLUMN = 0      # not a pivot, column from the block B
const UNKNOWN_PIVOT_COLUMN = 1  # maybe a pivot, undecided yet
const PIVOT_COLUMN = 2          # a pivot, column from the block A

mutable struct MacaulayMatrix{T <: Coeff}
    # A row in the matrix is represented (usually) sparsely by two vectors: the
    # support and the coefficients. The coefficients can be stored in different
    # places:
    # - Either as a pointer to a polynomial in the basis
    # - Or stored explicitly in the matrix, in this struct.

    upper_rows::Vector{Vector{ColumnLabel}}
    lower_rows::Vector{Vector{ColumnLabel}}

    # Explicitly stored coefficients of upper rows
    upper_coeffs::Vector{Vector{T}}
    # Explicitly stored coefficients of *some* rows. Usage of this field may
    # differ in runtime depending on the algorithm in use
    some_coeffs::Vector{Vector{T}}

    # Explicitly stored coefficients of block B.
    B_coeffs_dense::Vector{Vector{T}}
    # Explicitly stored coefficients of block D.
    D_coeffs_dense::Vector{Vector{T}}

    # Maps the columns of the matrix to corresponding monomial labels
    column_to_monom::Vector{MonomId}

    # If the AB part of the matrix is in Row-Reduced-Echelon-Form.
    # Generally, this field is `false`, and AB is only in Row-Echelon-Form.
    upper_part_is_rref::Bool

    # The number of columns in the left part of the matrix (in the blocks A, C)
    ncols_left::Int
    # The number of columns in the right part of the matrix (in the blocks B, D)
    ncols_right::Int

    # The number of rows actually filled in the upper/lower parts of the matrix
    nrows_filled_upper::Int
    nrows_filled_lower::Int

    # The number of pivots, i.e. the number of linearly independent rows
    # discovered in the block CD after performing linear reduction by AB
    npivots::Int
    # Maps a column to a pivot row. We have `pivots[i][1] == i`.
    # NOTE: `pivots[i]` may be un-assigned.
    pivots::Vector{Vector{ColumnLabel}}
    # A lightweight version of the above `pivots`, only for some specific cases.
    pivot_indices::Vector{Int}

    # Index of the row --> position of the array of coefficients for this row
    upper_to_coeffs::Vector{Int}
    lower_to_coeffs::Vector{Int}

    # Index of the row --> monomial multiplier
    upper_to_mult::Vector{MonomId}
    lower_to_mult::Vector{MonomId}

    # A vector of random elements from the ground field, which is used to
    # compute hashes of matrix rows
    hash_vector::Vector{T}

    sentinels::Vector{Int8}

    changematrix::Vector{Dict{Tuple{Int, MonomId}, T}}
end

# The number of allocated (not necessarily filled) rows and columns in the
# blocks in the matrix
function matrix_block_sizes(matrix::MacaulayMatrix)
    (
        length(matrix.upper_rows),
        length(matrix.lower_rows),
        matrix.ncols_left,
        matrix.ncols_right
    )
end

# The total number of allocated (not necessarily filled) rows and columns
function Base.size(matrix::MacaulayMatrix)
    mup, mlow, nleft, nright = matrix_block_sizes(matrix)
    (mup + mlow, nleft + nright)
end

# The number of actually filled rows in the matrix
function matrix_nrows_filled(matrix::MacaulayMatrix)
    (matrix.nrows_filled_upper, matrix.nrows_filled_lower)
end

# The number of actually filled columns in the matrix
function matrix_ncols_filled(matrix::MacaulayMatrix)
    (matrix.ncols_left, matrix.ncols_right)
end

# Initializes an empty matrix with coefficients of type T
function matrix_initialize(ring::PolyRing, ::Type{T}) where {T <: Coeff}
    MacaulayMatrix(
        Vector{Vector{ColumnLabel}}(),
        Vector{Vector{ColumnLabel}}(),
        Vector{Vector{T}}(),
        Vector{Vector{T}}(),
        Vector{Vector{T}}(),
        Vector{Vector{T}}(),
        Vector{MonomId}(),
        false,
        0,
        0,
        0,
        0,
        0,
        Vector{Vector{ColumnLabel}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{MonomId}(),
        Vector{MonomId}(),
        Vector{T}(),
        Vector{Int8}(),
        Vector{Dict{Tuple{Int, MonomId}, T}}()
    )
end

# Returns a string with human-readable information about the matrix. Currently,
# used in logging (call, e.g., `groebner(system, loglevel=-3)`)
function matrix_string_repr(matrix::MacaulayMatrix{T}) where {T}
    m, n = size(matrix)
    nupper, nlower = matrix_nrows_filled(matrix)
    nleft, nright = matrix_ncols_filled(matrix)
    m_A, n_A = nupper, nleft
    m_B, n_B = nupper, nright
    m_C, n_C = nlower, nleft
    m_D, n_D = nlower, nright
    nnz_A, nnz_B, nnz_C, nnz_D = 0, 0, 0, 0
    A_ref, A_rref = true, true
    # NOTE: probably want to adjust this when the terminal shrinks
    max_canvas_width = CustomUnicodePlots.DEFAULT_CANVAS_WIDTH
    canvas =
        CustomUnicodePlots.CanvasMatrix2x2(m_A, m_C, n_A, n_B, max_width=max_canvas_width)
    @inbounds for i in 1:nupper
        row = matrix.upper_rows[i]
        if row[1] < i
            A_ref = false
        end
        if row[1] != i
            A_rref = false
        end
        for j in 1:length(row)
            CustomUnicodePlots.point!(canvas, i, row[j])
            if row[j] <= nleft
                if row[j] != i
                    A_rref = false
                end
                nnz_A += 1
            else
                nnz_B += 1
            end
        end
    end
    @inbounds for i in 1:nlower
        row = matrix.lower_rows[i]
        for j in 1:length(row)
            CustomUnicodePlots.point!(canvas, nupper + i, row[j])
            if row[j] <= nleft
                nnz_C += 1
            else
                nnz_D += 1
            end
        end
    end
    nnz = nnz_A + nnz_B + nnz_C + nnz_D
    percent(x) = round(100 * x, digits=5)
    s = """
    $(typeof(matrix))
    $m x $n with $nnz nnz ($(percent(nnz / (m * n)))%)
    A: $(m_A) x $(n_A) with $(nnz_A) nnz ($(percent(nnz_A / (m_A * n_A)))%) (REF: $(A_ref), RREF: $(A_rref))
    B: $(m_B) x $(n_B) with $(nnz_B) nnz ($(percent(nnz_B / (m_B * n_B)))%)
    C: $(m_C) x $(n_C) with $(nnz_C) nnz ($(percent(nnz_C / (m_C * n_C)))%)
    D: $(m_D) x $(n_D) with $(nnz_D) nnz ($(percent(nnz_D / (m_D * n_D)))%)

    RREF flag: $(matrix.upper_part_is_rref)

    Sparsity pattern:

    $(canvas)
    """
    s
end

# It may be a good idea to call this on the entry to linear algebra
function matrix_well_formed(matrix::MacaulayMatrix)
    # TODO: not much is checked at the moment, but it can be improved :^)
    nupper, nlower = matrix_nrows_filled(matrix)
    nleft, nright = matrix_ncols_filled(matrix)
    true
end

function matrix_resize_upper_part!(matrix::MacaulayMatrix, size::Int)
    size <= length(matrix.upper_rows) && return nothing
    resize!(matrix.upper_rows, size)
    resize!(matrix.upper_to_coeffs, size)
    resize!(matrix.upper_to_mult, size)
    nothing
end

function matrix_resize_lower_part!(matrix::MacaulayMatrix, size::Int)
    size <= length(matrix.lower_rows) && return nothing
    resize!(matrix.lower_rows, size)
    resize!(matrix.lower_to_coeffs, size)
    resize!(matrix.lower_to_mult, size)
    nothing
end

# statistics_refresh and partially initialize the matrix
function matrix_reinitialize!(matrix::MacaulayMatrix, size::Int)
    new_size = size * 2
    matrix_resize_upper_part!(matrix, new_size)
    matrix_resize_lower_part!(matrix, new_size)
    matrix.upper_part_is_rref = false
    matrix.ncols_left = 0
    matrix.ncols_right = 0
    matrix.nrows_filled_upper = 0
    matrix.nrows_filled_lower = 0
    matrix.npivots = 0
    matrix
end

function matrix_resize_upper_part_if_needed!(matrix::MacaulayMatrix, new_size::Int)
    current_size = length(matrix.upper_rows)
    while current_size < new_size
        current_size *= 2
    end
    if current_size > length(matrix.upper_rows)
        matrix_resize_upper_part!(matrix, current_size)
    end
    nothing
end

###

###
# Change matrix

function matrix_changematrix_initialize!(matrix::MacaulayMatrix{T}, n) where {T <: Coeff}
    resize!(matrix.changematrix, n)
    for i in 1:length(matrix.changematrix)
        matrix.changematrix[i] = Dict{Tuple{Int, MonomId}, T}()
    end
    nothing
end

###
# Stuff

function matrix_convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix,
    basis::Basis{C},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters;
    batched_ht_insert::Bool=false
) where {C}
    # We mutate the basis directly by adding new elements
    @log :debug "# new pivots: $(matrix.npivots)"

    basis_resize_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.nprocessed

    _, _, nl, nr = matrix_block_sizes(matrix)
    support_size = sum(length, rows; init=0)
    @log :debug """
    After F4 reduction, in the right part of the matrix,
    # columns: $nr
    # pivots:  $(matrix.npivots)
    # support: $support_size
    """

    if batched_ht_insert && matrix.npivots > 1 && support_size > 4 * nr
        @log :debug "Using batched insert."
        column_to_basis_monom = zeros(MonomId, nr)
        @inbounds for i in 1:(matrix.npivots)
            for j in 1:length(rows[i])
                idx = rows[i][j]
                column_to_basis_monom[idx - nl] = 1
            end
        end
        matrix_insert_in_basis_hashtable_pivots_masked!(
            column_to_basis_monom,
            ht,
            symbol_ht,
            matrix.column_to_monom,
            nl
        )
        @inbounds for i in 1:(matrix.npivots)
            colidx = rows[i][1]
            for j in 1:length(rows[i])
                columnid = rows[i][j]
                rows[i][j] = column_to_basis_monom[columnid - nl]
            end
            basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
            basis.monoms[crs + i] = matrix.lower_rows[i]
            @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
        end
    else
        @inbounds for i in 1:(matrix.npivots)
            colidx = rows[i][1]
            matrix_insert_in_basis_hashtable_pivots!(
                rows[i],
                ht,
                symbol_ht,
                matrix.column_to_monom
            )
            basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
            basis.monoms[crs + i] = matrix.lower_rows[i]
            @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
        end
    end

    if params.changematrix
        resize!(basis.changematrix, basis.nfilled + matrix.npivots)
        for i in 1:(matrix.npivots)
            basis.changematrix[crs + i] = Dict{Int, Dict{MonomId, C}}()
            for ((poly_idx, poly_mult), cf) in matrix.changematrix[i]
                basis_changematrix_addmul!(
                    basis,
                    ht,
                    symbol_ht,
                    crs + i,
                    poly_idx,
                    poly_mult,
                    cf,
                    params.arithmetic
                )
            end
        end
    end

    basis.nfilled += matrix.npivots
end

function matrix_convert_rows_to_basis_elements_nf!(
    matrix::MacaulayMatrix,
    basis::Basis{C},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
) where {C}
    _, nlow = matrix_nrows_filled(matrix)
    basis_resize_if_needed!(basis, matrix.npivots)

    @inbounds for i in 1:nlow
        basis.nprocessed += 1
        basis.nnonredundant += 1
        basis.nonredundant[basis.nnonredundant] = basis.nprocessed
        if isassigned(matrix.some_coeffs, i)
            row = matrix.lower_rows[i]
            matrix_insert_in_basis_hashtable_pivots!(
                row,
                ht,
                symbol_ht,
                matrix.column_to_monom
            )
            basis.coeffs[basis.nprocessed] = matrix.some_coeffs[i]
            basis.monoms[basis.nprocessed] = row
        else
            empty!(basis.coeffs[basis.nprocessed])
            empty!(basis.monoms[basis.nprocessed])
        end
    end

    nothing
end

function matrix_polynomial_multiple_to_row!(
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable{M},
    basis_ht::MonomialHashtable{M},
    monom_hash::MonomHash,
    mult::M,
    poly::Vector{MonomId},
    skipfirst::Bool=false
) where {M <: Monom}
    row = similar(poly)
    hashtable_resize_if_needed!(symbol_ht, length(poly))

    hashtable_insert_polynomial_multiple!(
        row,
        monom_hash,
        mult,
        poly,
        basis_ht,
        symbol_ht,
        skipfirst
    )
end

function matrix_fill_column_to_monom_map!(
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable
)
    @invariant !symbol_ht.frozen

    hdata = symbol_ht.hashdata
    monoms = symbol_ht.monoms
    load = symbol_ht.load

    column_to_monom = Vector{MonomId}(undef, load - 1)
    j = 1
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        column_to_monom[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == PIVOT_COLUMN
            k += 1
        end
    end

    sort_partition_columns_by_labels!(column_to_monom, symbol_ht)

    cmp = let monoms = symbol_ht.monoms, ord = symbol_ht.ord
        function _cmp(x, y, ord)
            monom_isless(@inbounds(monoms[y]), @inbounds(monoms[x]), ord)
        end
        (x, y) -> _cmp(x, y, ord)
    end
    sort_part!(column_to_monom, 1, k, lt=cmp)
    sort_part!(column_to_monom, k + 1, length(column_to_monom), lt=cmp)

    matrix.ncols_left = k
    # -1 as long as hashtable load is always 1 more than actual
    matrix.ncols_right = load - matrix.ncols_left - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        hv = hdata[column_to_monom[k]]
        hdata[column_to_monom[k]] = Hashvalue(k, hv.hash, hv.divmask)
    end

    @inbounds for k in 1:(matrix.nrows_filled_upper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nrows_filled_lower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.column_to_monom = column_to_monom
end

# TODO: this does not belong here!...
function matrix_insert_in_basis_hashtable_pivots!(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomId}
) where {M <: Monom}
    hashtable_resize_if_needed!(ht, length(row))

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
        i = MonomHash(0)
        @inbounds while i <= ht.size
            k = hashtable_next_lookup_index(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if hashtable_is_hash_collision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        @invariant !ht.frozen

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask)

        ht.load += 1
    end

    nothing
end

# TODO: ...and this!
function matrix_insert_in_basis_hashtable_pivots_masked!(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomId},
    shift::Int
) where {M <: Monom}
    hashtable_resize_if_needed!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.monoms
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        if iszero(row[l])
            l += 1
            continue
        end

        hidx = column_to_monom[shift + l]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonomHash(0)
        @inbounds while i <= ht.size
            k = hashtable_next_lookup_index(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if hashtable_is_hash_collision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        @invariant !ht.frozen

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask)

        ht.load += 1
    end

    nothing
end
