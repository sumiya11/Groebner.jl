
#------------------------------------------------------------------------------

mutable struct Matrix{Tv}
    #=
        Matrix of the following structure

        | A  B |
        | C  D |

        A contains known pivots of reducing rows,
        and CD are rows to be reduced by AB

    =#

    # rows from upper, AB part of the matrix,
    # stored as vectors of corresponding exponents (already hashed)
    uprows::Vector{Vector{Int}}
    # rows from lower, CD part of the matrix,
    # stored as vectors of corresponding exponents (already hashed)
    lowrows::Vector{Vector{Int}}

    # maps column idx {1 ... ncols} to monomial hash {2 ... ht.load}
    # in some (?) hashtable
    col2hash::Vector{Int}

    # row coefficients
    # (some of the rows are stored in the basis,
    #  and some are stored here)
    coeffs::Vector{Vector{Tv}}

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

function initialize_matrix(ring::PolyRing{Tv}) where {Tv}
    uprows   = Vector{Vector{Int}}(undef, 0)
    lowrows  = Vector{Vector{Int}}(undef, 0)
    col2hash = Vector{Int}(undef, 0)
    coeffs   = Vector{Vector{Tv}}(undef, 0)

    size    = 0
    npivots = 0
    nrows   = 0
    ncols   = 0

    nup    = 0
    nlow   = 0
    nleft  = 0
    nright = 0

    up2coef   = Vector{Int}(undef, 0)
    low2coef   = Vector{Int}(undef, 0)

    Matrix(uprows, lowrows, col2hash, coeffs,
            size, npivots, nrows, ncols,
            nup, nlow, nleft, nright,
            up2coef, low2coef)
end

function reinitialize_matrix!(matrix, npairs)
    resize!(matrix.uprows, npairs*2)
    resize!(matrix.lowrows, npairs*2)
    resize!(matrix.up2coef, npairs*2)
    resize!(matrix.low2coef, npairs*2)
    matrix.size = 2*npairs
    matrix.ncols = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nup   = 0
    matrix.nlow  = 0
    matrix
end

#------------------------------------------------------------------------------

function normalize_sparse_row!(row)
    pinv = inv(row[1])
    for i in 2:length(row)
        row[i] *= pinv
    end
    row[1] = one(row[1])
    row
end

function reduce_dense_row_by_known_pivots_sparse!(
            dr, matrix, basis, pivs, dpiv, tmp_pos)

    ncols = matrix.ncols
    nleft = matrix.nleft
    k = 0
    np = -1

    for i in dpiv:ncols
        # if row element zero - no reduction
        if dr[i] == 0
            continue
        end
        # if pivot not defined
        if !isassigned(pivs, i)
            if np == -1
                np = i
            end
            k += 1
            continue
        end

        mul = inv(dr[i])
        dts = pivs[i]

        if i <= nleft
            cfs = basis.coeffs[matrix.up2coef[i]]
        else
            cfs = matrix.coeffs[matrix.low2coef[i]]
        end

        # TODO: Try doing better //
        for j in 1:length(dts)
            dr[dts[j]] -= mul * cfs[j]
        end
    end

    # all reduced !
    if k == 0
        return Vector{UInt}(undef, k)
    end

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    row = Vector{UInt}(undef, k)
    cf  = Vector{valtype(basis.coeffs[1])}(undef, k)
    j = 1
    for i in nleft:ncols
        if dr[i] != 0
            row[j] = i
            cf[j]  = dr[i]
            j += 1
        end
    end

    matrix.coeffs[tmp_pos] = cf

    return row
end

function exact_sparse_rref!(matrix, basis)
    ncols  = matrix.ncols
    nlow   = matrix.nlow
    nright = matrix.nright
    nleft  = matrix.nleft

    ground = parent(basis.coeffs[1][1])

    @info "entering sparse rref" matrix.nrows matrix.nlow matrix.nup
    @info "..."  ncols nright nleft
    @info ground

    # known pivots
    pivs = copy(matrix.uprows)

    # unknown pivots
    # (not discovered yet)
    upivs = matrix.lowrows

    col = 1

    for i in 1:nlow
        @debug "low row $i.."

        # select next row to be reduced
        npiv = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy)
        # we need to locate it somehow
        cfs  = basis.coeffs[matrix.up2coef[i]] # TODO

        k = 0

        # we load coefficients into dense array
        # TODO: move this
        drl = zeros(ground, ncols)
        for j in 1:length(npiv)
            drl[npiv[j]] = cfs[j]
        end

        # reduce it with known pivots from matrix.up
        npiv = reduce_dense_row_by_known_pivots_sparse!(drl, matrix, basis, pivs, col, i)
        # if fully reduced
        isempty(npiv) && continue

        if matrix.coeffs[i][1] != 1
            normalize_sparse_row!(matrix.coeffs[i])
        end
    end

    @debug "interreducing rows.."
    @info "matrix before::"
    dump(matrix, maxdepth=4)
    # number of new pivots
    npivs = 0

    # interreduce new pivots
    for i in 1:nright
        k = ncols - i
        @info "interrrrreducing " k
        if isassigned(matrix.lowrows, k)
            ds = matrix.lowrows[k]
            col = ds[1]
            dr = zeros(ground, ncols)
            cfs = matrix.coeffs[k]
            for j in 1:length(cfs)
                dr[ds[j]] = cfs[j]
            end
            npivs += 1
            matrix.lowrows[npivs] = reduce_dense_row_by_known_pivots_sparse!(dr, matrix, basis, pivs, col, 1)
        end
    end

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = npivs
    # resize!(matrix.lowrows, npivs)
end

function exact_sparse_linear_algebra!(matrix, basis)
    resize!(matrix.coeffs, matrix.nlow)

    @info "matrix in reduction"
    dump(matrix, maxdepth=3)

    exact_sparse_rref!(matrix, basis)
end


#------------------------------------------------------------------------------

function convert_hashes_to_columns!(matrix, symbol_ht)
    # col2hash = matrix.col2hash
    hdata    = symbol_ht.hashdata
    load     = symbol_ht.load

    @info "Converting hashes to columnds" symbol_ht.load
    # monoms from symbolic table represent one column in the matrix

    col2hash = Vector{UInt}(undef, load - 1)
    j = 1
    k = 0
    for i in symbol_ht.offset:load # TODO : hardcoding 2 is bad
        # map hash to column
        col2hash[j] = i
        j += 1

        # meaning the column is pivoted
        if hdata[i].idx == 2
            k += 1
        end
    end

    @info "init col2hash" col2hash

    # sort columns
    # TODO
    sort_columns_by_hash!(col2hash, symbol_ht)

    @info "after sort col2hash" col2hash

    matrix.nleft = k
    matrix.nright = load - matrix.nleft

    @info "updated matrix info" matrix.nleft matrix.nright

    # store the other direction of mapping,
    # hash -> column
    for k in 1:length(col2hash)
        hdata[col2hash[k]].idx = k
    end

    nterms = 0
    for k in 1:matrix.nup
        row = matrix.uprows[k]

        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end

        # TODO
        nterms += length(row)
    end

    for k in 1:matrix.nlow
        row = matrix.lowrows[k]

        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end

        # TODO
        nterms += length(row)
    end

    matrix.ncols = matrix.nleft + matrix.nright

    @assert matrix.nleft + matrix.nright == symbol_ht.load == matrix.ncols
    @assert matrix.nlow + matrix.nup == matrix.nrows

    matrix.col2hash = col2hash
end


function convert_matrix_rows_to_basis_elements!(
            matrix, basis, ht, symbol_ht)

    # we mutate basis array directly by adding new elements

    check_enlarge_basis!(basis, matrix.npivots)
    rows = matrix.lowrows
    crs = basis.ndone

    dump(matrix, maxdepth=4)

    for i in 1:matrix.npivots
        @info "converting matrix to basis" i rows[i]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.col2hash)
        # TODO : a constant

        basis.coeffs[crs + i] = matrix.coeffs[i]
        basis.gens[crs + i] = matrix.lowrows[i]
    end

    basis.ntotal += matrix.npivots
end
