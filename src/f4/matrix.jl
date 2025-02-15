# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Macaulay matrix

# MacaulayMatrix is a sparse matrix with columns labeled by distinct monomials.
# MacaulayMatrix represents a block matrix of the following structure:
#
#   | A  B |
#   | C  D |
#
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

    # Maps the columns of the matrix to corresponding monomial labels
    column_to_monom::Vector{MonomId}

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
    # May be un-assigned.
    pivots::Vector{Vector{ColumnLabel}}
    # A lightweight version of the above `pivots`, only for some specific cases.
    pivot_indices::Vector{Int}

    # Index of the row --> position of the array of coefficients for this row
    upper_to_coeffs::Vector{Int}
    lower_to_coeffs::Vector{Int}

    # Index of the row --> monomial multiplier
    upper_to_mult::Vector{MonomId}
    lower_to_mult::Vector{MonomId}

    sentinels::Vector{Int8}

    changematrix::Vector{Dict{Tuple{Int, MonomId}, T}}
end

# The number of allocated (not necessarily filled) rows and columns in the
# blocks in the matrix
function matrix_block_sizes(matrix::MacaulayMatrix)
    (matrix.nrows_filled_upper, matrix.nrows_filled_lower, matrix.ncols_left, matrix.ncols_right)
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
        Vector{MonomId}(),
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
        Vector{Int8}(),
        Vector{Dict{Tuple{Int, MonomId}, T}}()
    )
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

function matrix_reinitialize!(matrix::MacaulayMatrix, size::Int)
    new_size = size * 2
    matrix_resize_upper_part!(matrix, new_size)
    matrix_resize_lower_part!(matrix, new_size)
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

    basis_resize_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.n_processed

    _, _, nl, nr = matrix_block_sizes(matrix)

    @invariant all(i -> isdefined(rows, i), 1:length(rows))
    support_size = sum(length, rows; init=0)

    if batched_ht_insert && matrix.npivots > 1 && support_size > 4 * nr
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
            matrix_insert_in_basis_hashtable_pivots!(rows[i], ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
            basis.monoms[crs + i] = matrix.lower_rows[i]
            @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
        end
    end

    if params.changematrix
        resize!(basis.changematrix, basis.n_filled + matrix.npivots)
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

    basis.n_filled += matrix.npivots
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
        basis.n_processed += 1
        basis.n_nonredundant += 1
        basis.nonredundant_indices[basis.n_nonredundant] = basis.n_processed
        if isassigned(matrix.some_coeffs, i)
            row = matrix.lower_rows[i]
            matrix_insert_in_basis_hashtable_pivots!(row, ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[basis.n_processed] = matrix.some_coeffs[i]
            basis.monoms[basis.n_processed] = row
        else
            empty!(basis.coeffs[basis.n_processed])
            empty!(basis.monoms[basis.n_processed])
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
    poly::Vector{MonomId}
) where {M <: Monom}
    row = similar(poly)
    hashtable_resize_if_needed!(symbol_ht, length(poly))

    hashtable_insert_polynomial_multiple!(row, monom_hash, mult, poly, basis_ht, symbol_ht)
end

function matrix_fill_column_to_monom_map!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)
    @invariant !symbol_ht.frozen

    column_to_monom = Vector{MonomId}(undef, symbol_ht.load - 1)
    j = 1
    k = 0
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        column_to_monom[j] = i
        j += 1
        if symbol_ht.labels[i] == PIVOT_COLUMN
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
    # -1 as long as hashtable symbol_ht.load is always 1 more than actual
    matrix.ncols_right = symbol_ht.load - matrix.ncols_left - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        symbol_ht.labels[column_to_monom[k]] = k
    end

    @inbounds for k in 1:(matrix.nrows_filled_upper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = symbol_ht.labels[row[j]]
        end
    end

    @inbounds for k in 1:(matrix.nrows_filled_lower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = symbol_ht.labels[row[j]]
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

    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
    bexps = ht.monoms
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = column_to_monom[row[l]]
        h = symbol_ht.hashvals[hidx]

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

        ht.labels[pos] = symbol_ht.labels[hidx]
        ht.hashvals[pos] = h
        ht.divmasks[pos] = symbol_ht.divmasks[hidx]

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

    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
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

        h = symbol_ht.hashvals[hidx]

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

        ht.labels[pos] = symbol_ht.labels[hidx]
        ht.hashvals[pos] = h
        ht.divmasks[pos] = symbol_ht.divmasks[hidx]

        ht.load += 1
    end

    nothing
end
