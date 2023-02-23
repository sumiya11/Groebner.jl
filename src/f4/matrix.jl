# Linear algebra and MacaulayMatrix

struct LinearAlgebraContent

mutable struct MacaulayMatrix{T<:Coeff}
    #=
        Matrix of the structure

        | A  B |
        | C  D |

        part A contains known pivots of reducing rows,
        and part CD contains rows to be reduced by AB
    =#

    # rows from upper, AB part of the matrix,
    # stored as vectors of corresponding exponents (already hashed and sorted)
    uprows::Vector{Vector{ColumnIdx}}
    # rows from lower, CD part of the matrix,
    # stored as vectors of corresponding exponents (already hashed and sorted)
    lowrows::Vector{Vector{ColumnIdx}}

    # maps column idx {1 ... ncols} to monomial hash idx {2 ... ht.load}
    # in some hashtable
    col2hash::Vector{MonomIdx}

    # row coefficients
    # (some of the rows are stored in the basis,
    # and some are stored here)
    coeffs::Vector{Vector{T}}

    #= sizes info =#
    # total number of allocated rows
    size::Int
    # number of pivots,
    # ie new basis elements discovered after matrix reduction
    npivots::Int
    # number of filled rows, nrows <= size
    nrows::Int
    # number of columns
    ncols::Int

    # number of upper rows (in AB section)
    nup::Int
    # number of lower rows (in CD section)
    nlow::Int
    # number of left cols  (in AC section)
    nleft::Int
    # number of right cols (in BD section)
    nright::Int

    # maps column idx {1 ... ncols} to index of coefficient array
    # First nleft indices point to coefficients from AB part of the matrix
    # (these coefficients are owned by the basis struct)
    # Last nright indices point to coefficients from CD part of the matrix
    # (these are owned by this object and stored in coeffs array)
    # Essentially each coefficient array from the first part represents a reducer row,
    # and each array from the second part stands for a reduced *nonzero* row
    # should be row to coef
    up2coef::Vector{Int}
    low2coef::Vector{Int}
end

# Initializes an empty matrix with coefficients of type T
function initialize_matrix(ring::PolyRing, ::Type{T}) where {T<:Coeff}
    uprows = Vector{Vector{ColumnIdx}}(undef, 0)
    lowrows = Vector{Vector{ColumnIdx}}(undef, 0)
    col2hash = Vector{MonomIdx}(undef, 0)
    coeffs = Vector{Vector{T}}(undef, 0)

    size = 0
    npivots = 0
    nrows = 0
    ncols = 0

    nup = 0
    nlow = 0
    nleft = 0
    nright = 0

    up2coef = Vector{Int}(undef, 0)
    low2coef = Vector{Int}(undef, 0)

    MacaulayMatrix(uprows, lowrows, col2hash, coeffs,
        size, npivots, nrows, ncols,
        nup, nlow, nleft, nright,
        up2coef, low2coef)
end

# Refresh and initialize matrix for `npairs` elements
function reinitialize_matrix!(matrix::MacaulayMatrix{T}, npairs::Int) where {T}
    resize!(matrix.uprows, npairs * 2)
    resize!(matrix.lowrows, npairs * 2)
    resize!(matrix.up2coef, npairs * 2)
    resize!(matrix.low2coef, npairs * 2)
    matrix.size = 2 * npairs
    matrix.ncols = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nup = 0
    matrix.nlow = 0
    matrix
end

#------------------------------------------------------------------------------

function linear_algebra!(ring::PolyRing, matrix::MacaulayMatrix, basis::Basis, linalg::Val{:exact}, rng)
    resize!(matrix.coeffs, matrix.nlow)
    exact_sparse_rref!(ring, matrix, basis)
end

function linear_algebra!(ring::PolyRing, matrix::MacaulayMatrix, basis::Basis, linalg::Val{:prob}, rng)
    resize!(matrix.coeffs, matrix.nlow)
    randomized_sparse_rref!(ring, matrix, basis, rng)
end

#------------------------------------------------------------------------------

# if the field is integers modulo prime, enable optimization
function select_divisor(coeffs::Vector{Vector{T}}, ch) where {T<:CoeffFF}
    Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(ch)
end

function select_divisor(coeffs::Vector{Vector{T}}, ch) where {T<:CoeffQQ}
    ch
end

# Normalize `row` by first coefficient
#
# Finite field magic specialization
function normalize_sparse_row!(row::Vector{T}, magic) where {T<:CoeffFF}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = invmod(row[1], magic.divisor) % magic
    @inbounds for i in 2:length(row)
        row[i] = (row[i] * pinv) % magic
    end
    @inbounds row[1] = one(row[1])
    row
end

# Normalize `row` by first coefficient
#
# Rational field specialization
function normalize_sparse_row!(row::Vector{T}, ch) where {T<:CoeffQQ}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = inv(row[1])
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end
    @inbounds row[1] = one(row[1])
    row
end

function ui_redmod_barrett_half!(row, indices, cfs, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    p = magic.divisor
    ϕ = div(mul << 32, p, RoundUp)
    @inbounds for i in 1:length(indices)
        idx = indices[i]
        c = (cfs[i] * ϕ) >> 32
        x = row[idx] + mul*cfs[i] - c*p
        row[idx] = min(x, x - p)
    end
    row
end

# This needs further tuning, not used at the moment
function u64_red_inline_barrett_half_4t32!(row::Vector{T}, indices, cfs, magic) where {T<:UInt32}    
    corr = isodd(length(indices))
    lastidx = length(indices) - corr

    N = 4
    @inbounds mul = magic.divisor - row[indices[1]]
    mulv = SIMD.Vec{N, UInt32}(mul)
    pv = SIMD.Vec{N, UInt32}(magic.divisor)
    ϕv = SIMD.Vec{2, UInt64}(div(UInt64(mul) << 32, magic.divisor))
    zerov = SIMD.Vec{2, UInt64}(0)

    @inbounds for i in 1:N:lastidx
        cfs1, cfs2, cfs3, cfs4 = cfs[i], cfs[i + 1], cfs[i + 2], cfs[i + 3]
        cfs12 = SIMD.Vec{2, UInt64}((cfs1, cfs2))
        cfs34 = SIMD.Vec{2, UInt64}((cfs3, cfs4))
        
        c12 = cfs12 * ϕv
        c34 = cfs34 * ϕv
        c12 = unpack_hi(c12, zerov)
        c34 = unpack_hi(c34, zerov)

        c1234 = SIMD.Vec{N, UInt32}((c12[1], c12[2], c34[1], c34[2]))

        cfs1234 = SIMD.Vec{4, UInt32}((cfs1, cfs2, cfs3, cfs4))
        idx1234 = vload(SIMD.Vec{N, Int}, indices, i)
        row1234 = vgather(row, idx1234)
        
        x1234 = row1234 + cfs1234 * mulv - c1234*pv
        x1234 = min(x1234, x1234 - pv)
        
        SIMD.vscatter(x1234, row, idx1234)
    end

    row
end

# reduces row by mul*cfs modulo ch at indices positions
#
# Finite field magic specialization
function reduce_by_pivot!(row::Vector{T}, indices::Vector{ColumnIdx},
        cfs::Vector{T}, magic::Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{T}) where {T<:CoeffFF}

    @inbounds mul = magic.divisor - row[indices[1]]

    # on our benchmarks usually
    # length(row) / length(indices) varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + mul * cfs[j]) % magic
    end

    nothing
end

# reduces row by mul*cfs at indices positions
#
# Rational field specialization
function reduce_by_pivot!(row::Vector{T}, indices::Vector{ColumnIdx},
        cfs::Vector{T}, ch) where {T<:CoeffQQ}

    @inbounds mul = -row[indices[1]]

    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * cfs[j]
    end

    nothing
end

# zero entries of densecoeffs and load coefficients cfsref to indices rowexps
#
# Finite field specialization
function load_indexed_coefficients!(densecoeffs::Vector{T}, rowexps, cfsref) where {T<:CoeffFF}
    @inbounds for i in 1:length(densecoeffs)
        densecoeffs[i] = T(0)
    end
    @inbounds for j in 1:length(rowexps)
        densecoeffs[rowexps[j]] = cfsref[j]
    end
end

# zero entries of densecoeffs and load coefficients cfsref to indices rowexps
#
# Rational numbers specialization
function load_indexed_coefficients!(densecoeffs::Vector{T}, rowexps, cfsref) where {T<:CoeffQQ}
    densecoeffs .= T(0)
    @inbounds for j in 1:length(rowexps)
        densecoeffs[rowexps[j]] = cfsref[j]
    end
end

#------------------------------------------------------------------------------

function reduce_dense_row_by_known_pivots_sparse!(
        densecoeffs::Vector{C}, matrix::MacaulayMatrix{C}, basis::Basis{C},
        pivs::Vector{Vector{ColumnIdx}}, startcol, tmp_pos, magic;
        exact_colmap::Bool=false) where {C<:Coeff}
    reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, 
            ColumnIdx(startcol), ColumnIdx(tmp_pos), magic; 
            exact_colmap=exact_colmap)
end

const counter = Ref{Int}(0)

function reduce_dense_row_by_known_pivots_sparse!(
    densecoeffs::Vector{C}, matrix::MacaulayMatrix{C}, basis::Basis{C},
    pivs::Vector{Vector{ColumnIdx}}, startcol::ColumnIdx, tmp_pos::ColumnIdx, magic;
    exact_colmap::Bool=false) where {C<:Coeff}

    ncols = matrix.ncols
    nleft = matrix.nleft

    # new row nonzero elements count
    k = 0
    uzero = zero(C)

    # new pivot index
    np = -1

    counter[] = 0

    @inbounds for i in startcol:ncols
        # if row element zero - no reduction
        if iszero(densecoeffs[i])
            continue
        end

        if !isassigned(pivs, i) || (tmp_pos != -1 && tmp_pos == i)
            if np == -1
                np = i
            end
            k += 1
            continue
        end

        # exponents of reducer row at column i
        reducerexps = pivs[i]

        if exact_colmap # if reducer from new matrix pivots
            @inbounds cfs = matrix.coeffs[tmp_pos]
        elseif i <= nleft # if reducer is from upper part of the matrix
            @inbounds cfs = basis.coeffs[matrix.up2coef[i]]
        else # if reducer is from upper part of the matrix
            @inbounds cfs = matrix.coeffs[matrix.low2coef[i]]
        end
        counter[] += 1

        reduce_by_pivot!(densecoeffs, reducerexps, cfs, magic)
    end
    println(length(densecoeffs), ", ", ", ", counter[])

    newrow = Vector{ColumnIdx}(undef, k)
    newcfs = Vector{C}(undef, k)

    # all reduced !
    if k == 0
        return true, newrow, newcfs
    end

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    @assert k > 0
    j = 1
    @inbounds for i in np:ncols # starting from new pivot
        if densecoeffs[i] != uzero
            newrow[j] = i
            newcfs[j] = densecoeffs[i]
            j += 1
        end
    end

    return false, newrow, newcfs
end

#------------------------------------------------------------------------------

function absolute_pivots!(matrix::MacaulayMatrix)
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:matrix.nup
        pivs[i] = matrix.uprows[i]
    end
    l2c_tmp = Vector{ColumnIdx}(undef, max(matrix.ncols, matrix.nlow))
    @inbounds for i in 1:matrix.nlow
        l2c_tmp[matrix.lowrows[i][1]] = matrix.low2coef[i]
    end
    pivs, l2c_tmp
end

function interreduce_lower_part!(matrix::MacaulayMatrix{C}, basis::Basis{C}, pivs, magic; reversed=false) where {C<:Coeff}
    # number of new pivots
    newpivs = 0
    densecoeffs = zeros(C, matrix.ncols)
    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)
    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:matrix.nright
        k = matrix.ncols - i + 1
        !isassigned(pivs, k) && continue
        if k <= matrix.nleft # upper part of matrix
            cfsref = basis.coeffs[matrix.up2coef[k]]
        else # lower part of matrix
            cfsref = matrix.coeffs[matrix.low2coef[k]]
        end

        @assert length(cfsref) == length(pivs[k])

        startcol = pivs[k][1]
        load_indexed_coefficients!(densecoeffs, pivs[k], cfsref)

        _, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, startcol, magic)
        newpivs += 1

        # update row and coeffs
        if !reversed
            matrix.lowrows[newpivs] = newrow
            matrix.coeffs[matrix.low2coef[k]] = newcfs
            pivs[k] = matrix.lowrows[newpivs]
        else
            matrix.lowrows[matrix.nrows - newpivs + 1] = newrow
            matrix.coeffs[matrix.low2coef[k]] = newcfs
            pivs[k] = matrix.lowrows[matrix.nrows - newpivs + 1]
        end
    end

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = newpivs
    resize!(matrix.lowrows, newpivs)
end

# Linear algebra option 1
function exact_sparse_rref!(ring, matrix::MacaulayMatrix{C}, 
                        basis::Basis{C}) where {C<:Coeff}
    ncols = matrix.ncols
    nlow = matrix.nlow

    magic = select_divisor(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    # YES
    pivs, l2c_tmp = absolute_pivots!(matrix)
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

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
        load_indexed_coefficients!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1, magic)
        # if fully reduced
        zeroed && continue

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.low2coef[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], magic)
    end

    interreduce_lower_part!(matrix, basis, pivs, magic)
end

function nblocks_in_randomized(nrows::Int)
    floor(Int, sqrt(nrows / 3)) + 1
end

function randomized_sparse_rref!(ring::PolyRing, matrix::MacaulayMatrix{C}, 
                        basis::Basis{C}, rng) where {C<:Coeff}
    ncols = matrix.ncols
    nlow = matrix.nlow

    magic = select_divisor(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = absolute_pivots!(matrix)
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

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
                mulcoeffs[j] = rand(rng, C) % magic
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
                    densecoeffs[ridx] = (densecoeffs[ridx] + mulcoeffs[k] * cfsref[l]) % magic
                end
            end

            # reduce it with known pivots from matrix.uprows
            # first nonzero in densecoeffs is at startcol position
            zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1, magic)

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
            matrix.low2coef[newrow[1]] = absolute_i

            # normalize if needed
            normalize_sparse_row!(matrix.coeffs[absolute_i], magic)
            
            ctr += 1
        end
    end

    interreduce_lower_part!(matrix, basis, pivs, magic)
end

function exact_sparse_rref_interreduce!(ring::PolyRing, matrix::MacaulayMatrix{C}, basis::Basis{C}) where {C<:Coeff}
    resize!(matrix.lowrows, matrix.ncols)
    resize!(matrix.up2coef, matrix.ncols)
    resize!(matrix.low2coef, matrix.ncols)
    resize!(matrix.coeffs, matrix.ncols)

    magic = select_divisor(matrix.coeffs, ring.ch)

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:matrix.nrows
        pivs[matrix.uprows[i][1]] = matrix.uprows[i]
        matrix.low2coef[matrix.uprows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.up2coef[i]])
    end

    interreduce_lower_part!(matrix, basis, pivs, magic, reversed=true)
end

function exact_sparse_rref_isgroebner!(ring::PolyRing, matrix::MacaulayMatrix{C}, basis::Basis{C}) where {C<:Coeff}
    magic = select_divisor(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = absolute_pivots!(matrix)
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    upivs = matrix.lowrows

    densecoeffs = zeros(C, matrix.ncols)

    @inbounds for i in 1:matrix.nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]
        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]
        load_indexed_coefficients!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, _, _ = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1, magic)
        # @warn "reduced " zeroed newrow newcfs
        # if fully reduced
        zeroed && continue
        return false
    end
    return true
end

function exact_sparse_rref_nf!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    tobereduced::Basis{C},
    basis::Basis{C}) where {C<:Coeff}

    resize!(matrix.coeffs, matrix.nlow)
    
    ncols = matrix.ncols
    nlow = matrix.nlow

    magic = select_divisor(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = absolute_pivots!(matrix)
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    upivs = matrix.lowrows

    densecoeffs = zeros(C, ncols)

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = tobereduced.coeffs[rowidx2coef[i]]
        load_indexed_coefficients!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        # zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1, magic)
        # @warn "reduced " zeroed
        # if fully reduced
        zeroed && continue

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        matrix.lowrows[i] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.low2coef[i] = i
    end
    matrix.npivots = matrix.nrows = matrix.size = matrix.nlow
end

#------------------------------------------------------------------------------

function column_to_monom_mapping!(
        matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    col2hash = Vector{MonomIdx}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    @inbounds for i in symbol_ht.offset:load
        # column to hash index
        col2hash[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == 2
            k += 1
        end
    end

    sort_columns_by_hash!(col2hash, symbol_ht)

    matrix.nleft = k
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(col2hash)
        hdata[col2hash[k]].idx = k
    end

    @inbounds for k in 1:matrix.nup
        row = matrix.uprows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:matrix.nlow
        row = matrix.lowrows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.ncols = matrix.nleft + matrix.nright

    @assert matrix.nleft + matrix.nright == symbol_ht.load - 1 == matrix.ncols
    @assert matrix.nlow + matrix.nup == matrix.nrows

    matrix.col2hash = col2hash
end

function convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix, basis::Basis,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # we mutate basis array directly by adding new elements

    check_enlarge_basis!(basis, matrix.npivots)
    rows = matrix.lowrows
    crs = basis.ndone

    @inbounds for i in 1:matrix.npivots
        colidx = rows[i][1]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.col2hash)
        basis.coeffs[crs+i] = matrix.coeffs[matrix.low2coef[colidx]]
        basis.monoms[crs+i] = matrix.lowrows[i]
    end

    basis.ntotal += matrix.npivots
end

function convert_rows_to_basis_elements_nf!(
    matrix::MacaulayMatrix, basis::Basis,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    check_enlarge_basis!(basis, matrix.npivots)

    @inbounds for i in 1:matrix.npivots
        basis.ndone += 1
        basis.nlead += 1
        basis.nonred[basis.nlead] = basis.ndone
        if isassigned(matrix.coeffs, i)
            row = matrix.lowrows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.col2hash)
            basis.coeffs[basis.ndone] = matrix.coeffs[i]
            basis.monoms[basis.ndone] = row
        else
            empty!(basis.coeffs[basis.ndone])
            empty!(basis.monoms[basis.ndone])
        end
    end

    nothing
end
