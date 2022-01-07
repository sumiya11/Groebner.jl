
#------------------------------------------------------------------------------

mutable struct MacaulayMatrix
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
    coeffs::Vector{Vector{UInt64}}

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

function initialize_matrix(ring::PolyRing)
    uprows   = Vector{Vector{Int}}(undef, 0)
    lowrows  = Vector{Vector{Int}}(undef, 0)
    col2hash = Vector{Int}(undef, 0)
    coeffs   = Vector{Vector{UInt64}}(undef, 0)

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

    MacaulayMatrix(uprows, lowrows, col2hash, coeffs,
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

function normalize_sparse_row!(row::Vector{UInt64}, ch::UInt64)
    pinv = uinvmod(row[1], ch)
    for i in 2:length(row)
        # row[i] *= pinv
        row[i] = umultmod(row[i], pinv, ch)
    end
    row[1] = one(row[1])
    row
end

function reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs::Vector{UInt64}, matrix::MacaulayMatrix, basis::Basis,
            pivs, startcol, tmp_pos)

    @info "reducing $densecoeffs by"
    #println("pivs = ", pivs)

    ncols = matrix.ncols
    nleft = matrix.nleft

    ch = basis.ch

    # new row nonzero elements count
    k = 0

    # new pivot index
    np = -1

    for i in startcol:ncols
        # if row element zero - no reduction
        if densecoeffs[i] == 0
            continue
        end
        # if pivot not defined
        #= WARNING =#

        if !isassigned(pivs, i) || (tmp_pos != -1 && tmp_pos == i)
            if np == -1
                np = i
            end
            k += 1
            continue
        end


        # mul = -densecoeffs[i]
        mul = ucompmod(densecoeffs[i], ch)

        # exponents of reducer row at column i
        reducerexps = pivs[i]

        # if reducer is from upper
        if i <= nleft
            cfs = basis.coeffs[matrix.up2coef[i]]
        else # of lower part of matrix
            cfs = matrix.coeffs[matrix.low2coef[i]]
        end

        #println(reducerexps)
        #println(cfs)

        for j in 1:length(reducerexps)
            # densecoeffs[reducerexps[j]] += mul * cfs[j]
            densecoeffs[reducerexps[j]] = umultsummod(densecoeffs[reducerexps[j]], mul, cfs[j], ch)
        end
    end

    newrow = Vector{Int}(undef, k)
    newcfs = Vector{UInt64}(undef, k)

    # all reduced !
    if k == 0
        return true, newrow, newcfs
    end

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    j = 1
    for i in np:ncols # from new pivot
        if densecoeffs[i] != 0
            newrow[j] = i
            newcfs[j] = densecoeffs[i]
            j += 1
        end
    end

    return false, newrow, newcfs
end

function exact_sparse_rref!(matrix::MacaulayMatrix, basis::Basis)
    ncols  = matrix.ncols
    nlow   = matrix.nlow
    nright = matrix.nright
    nleft  = matrix.nleft

    # TODO
    # ground = parent(basis.coeffs[1][1])

    @info "entering sparse rref" matrix.nrows matrix.nlow matrix.nup
    @info "..."  ncols nright nleft
    
    # known pivots
    pivs = Vector{Vector{Int}}(undef, ncols)
    for i in 1:matrix.nup
        pivs[i] = copy(matrix.uprows[i])
    end

    #=
    WARNing. TODO
    =#
    #=
    pivs = Vector{Vector{Int}}(undef, ncols)
    up2coef = Vector{Int}(undef, ncols)
    for ui in 1:matrix.nup
        pivs[matrix.uprows[ui][1]] = copy(matrix.uprows[ui])
        up2coef[matrix.uprows[ui][1]] = matrix.up2coef[ui]
    end
    matrix.up2coef = up2coef
    =#
    #=
    WARNing. TODO
    =#

    resize!(matrix.low2coef, ncols)

    @info "known pivots" pivs

    # unknown pivots
    # (not discovered yet)
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

    @info "unknown pivots" upivs

    rowcoeffs = zeros(UInt64, ncols)
    for i in 1:nlow
        @debug "low row $i.."

        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref  = basis.coeffs[matrix.low2coef[i]]

        k = 0

        # we load coefficients into dense array
        # into rowexps indices
        # TODO: move this
        densecoeffs = zeros(UInt64, ncols)
        for j in 1:length(rowexps)
            densecoeffs[rowexps[j]] = cfsref[j]
        end

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        @debug "low row $i reduced into" zeroed newrow newcfs

        # if fully reduced
        zeroed && continue

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        if !isassigned(pivs, newrow[1])
            pivs[newrow[1]]  = newrow
            # set ref to coefficient to matrix
            # guaranteed to be from lower part
            matrix.low2coef[newrow[1]] = i
        end

        # normalize if needed
        if matrix.coeffs[i][1] != 1
            normalize_sparse_row!(matrix.coeffs[i], basis.ch)
        end
    end

    @debug "interreducing rows.."
    @info "matrix before::"
    #dump(matrix, maxdepth=4)
    #println("PIVS ARE : \n", pivs)
    # number of new pivots
    newpivs = 0

    densecfs = zeros(UInt64, ncols)
    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)

    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    for i in 1:nright
        k = ncols - i + 1
        @info "interrrrreducing right column" k
        if isassigned(pivs, k)
            @info "pivot $k exists"

            densecfs = zeros(UInt64, ncols)

            if k <= nleft
                cfsref = basis.coeffs[matrix.up2coef[k]]
            else # of lower part of matrix
                cfsref = matrix.coeffs[matrix.low2coef[k]]
            end

            startcol = pivs[k][1]

            @assert length(cfsref) == length(pivs[k])

            for j in 1:length(pivs[k])
                densecfs[pivs[k][j]] = cfsref[j]
            end
            newpivs += 1

            @debug "reducing" densecfs startcol

            zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecfs, matrix, basis, pivs, startcol, startcol)
            @debug "reduced" zeroed newrow newcfs
            # TODO: fails
            # @assert !zeroed

            # update row and coeffs
            matrix.lowrows[newpivs] = newrow
            matrix.coeffs[matrix.low2coef[k]] = newcfs
            pivs[k] = matrix.lowrows[newpivs]
        end
    end

    @debug "INTERREDUCED" newpivs
    @info "matrix is::"
    #dump(matrix, maxdepth=4)
    #println("PIVS ARE : \n", pivs)

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = newpivs
    resize!(matrix.lowrows, newpivs)
end

function exact_sparse_linear_algebra!(matrix::MacaulayMatrix, basis::Basis)
    resize!(matrix.coeffs, matrix.nlow)

    @info "matrix in reduction"
    #dump(matrix, maxdepth=3)

    exact_sparse_rref!(matrix, basis)
end


#------------------------------------------------------------------------------

function convert_hashes_to_columns!(
            matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # col2hash = matrix.col2hash
    hdata    = symbol_ht.hashdata
    load     = symbol_ht.load

    @info "Converting hashes to columnds" symbol_ht.load
    # monoms from symbolic table represent one column in the matrix

    col2hash = Vector{UInt}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    for i in symbol_ht.offset:load
        # column to hash index
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

    @info "hashtable"
    #println(symbol_ht.exponents, "\nload: ", symbol_ht.load)
    #println(symbol_ht.hashdata)

    @info "after sort col2hash" col2hash

    matrix.nleft = k
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

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

    @assert matrix.nleft + matrix.nright == symbol_ht.load-1 == matrix.ncols
    @assert matrix.nlow + matrix.nup == matrix.nrows

    matrix.col2hash = col2hash
end


function convert_matrix_rows_to_basis_elements!(
            matrix, basis, ht, symbol_ht)

    # we mutate basis array directly by adding new elements

    check_enlarge_basis!(basis, matrix.npivots)
    rows = matrix.lowrows
    crs = basis.ndone

    #dump(matrix, maxdepth=4)

    @info "before insert_in_basis_hash_table_pivots"
    #println(ht.exponents, "\n load: ", ht.load)
    @info "symbolic::"
    #println(symbol_ht.exponents, "\n symbol_ht: ", symbol_ht.load)

    @info "col2hash" matrix.col2hash

    for i in 1:matrix.npivots
        @info "converting matrix to basis" i rows[i]
        colidx = rows[i][1]

        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.col2hash)
        # TODO : a constant

        # an interesing way to find coefficients
        basis.coeffs[crs + i] = matrix.coeffs[matrix.low2coef[colidx]]
        basis.gens[crs + i] = matrix.lowrows[i]
        # sort_gens_terms_decreasing!(basis, ht, crs + i)
    end

    basis.ntotal += matrix.npivots
end
