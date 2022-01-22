
#------------------------------------------------------------------------------

function normalize_sparse_row!(row::Vector{UInt64}, ch::UInt64)
    pinv = uinvmod(row[1], ch)
    @inbounds for i in 2:length(row)
        # row[i] *= pinv
        row[i] = umultmod(row[i], pinv, ch)
    end
    row[1] = one(row[1])
    row
end

# reduces row by mul*cfs at indices positions
function reduce_by_pivot!(row::Vector{UInt64}, indices::Vector{Int},
            mul::UInt64, cfs::Vector{UInt64}, ch::UInt64)

    # length(row) / length(indices) varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + mul*cfs[j]) % ch
    end
    nothing
end

function reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs::Vector{UInt64}, matrix::MacaulayMatrix, basis::Basis,
            pivs::Vector{Vector{Int}}, startcol::Int, tmp_pos::Int)

    ncols = matrix.ncols
    nleft = matrix.nleft

    ch = basis.ch

    # new row nonzero elements count
    k = 0
    uzero = UInt64(0)

    # new pivot index
    np = -1

    for i in startcol:ncols
        # if row element zero - no reduction
        @inbounds if densecoeffs[i] == uzero
            continue
        end
        # if pivot not defined
        #= WARNING =#

        # @warn "inner iter $i" startcol tmp_pos pivs
        if !isassigned(pivs, i) || (tmp_pos != -1 && tmp_pos == i)
            # println("skipped")
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
            @inbounds cfs = basis.coeffs[matrix.up2coef[i]]
        else # of lower part of matrix
            @inbounds cfs = matrix.coeffs[matrix.low2coef[i]]
        end

        reduce_by_pivot!(densecoeffs, reducerexps, mul, cfs, ch)
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
    @inbounds for i in np:ncols # from new pivot
        @inbounds if densecoeffs[i] != uzero
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

    # known pivots
    # no_copy
    pivs = Vector{Vector{Int}}(undef, ncols)
    @inbounds for i in 1:matrix.nup
        pivs[i] = matrix.uprows[i]
    end

    # CHANGED in order to prevent bug
    # when several rows in the matrix are equal
    l2c_tmp = Vector{Int}(undef, max(ncols, matrix.nlow))
    @inbounds for i in 1:nlow
        l2c_tmp[matrix.lowrows[i][1]] = matrix.low2coef[i]
    end
    # CHANGED
    # no_copy
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots
    # (not discovered yet)
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

    densecoeffs = zeros(UInt64, ncols)

    #=
    @warn "before reducing low"
    dump(matrix, maxdepth=5)
    @warn "lowrow2coef" rowidx2coef
    =#

    for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref  = basis.coeffs[rowidx2coef[i]]

        k = 0

        # we load coefficients into dense array
        # into rowexps indices
        # TODO: move this
        densecoeffs .= UInt64(0)
        @inbounds for j in 1:length(rowexps)
            densecoeffs[rowexps[j]] = cfsref[j]
        end

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        # zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        # @warn "reduced " zeroed newrow newcfs
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

    #=
    @warn "after reducing low"
    dump(matrix, maxdepth=5)
    println("PIVS: ", pivs)
    =#

    # number of new pivots
    newpivs = 0

    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)

    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    densecfs = zeros(UInt64, ncols)

    for i in 1:nright
        k = ncols - i + 1
        if isassigned(pivs, k)
            densecfs .= UInt64(0)

            if k <= nleft
                cfsref = basis.coeffs[matrix.up2coef[k]]
            else # of lower part of matrix
                cfsref = matrix.coeffs[matrix.low2coef[k]]
            end

            startcol = pivs[k][1]

            @assert length(cfsref) == length(pivs[k])

            @inbounds for j in 1:length(pivs[k])
                densecfs[pivs[k][j]] = cfsref[j]
            end
            newpivs += 1

            zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecfs, matrix, basis, pivs, startcol, startcol)

            # update row and coeffs
            matrix.lowrows[newpivs] = newrow
            matrix.coeffs[matrix.low2coef[k]] = newcfs
            matrix.low2coef[k] = matrix.low2coef[k]
            pivs[k] = matrix.lowrows[newpivs]
        end
    end

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = newpivs
    resize!(matrix.lowrows, newpivs)
end

function exact_sparse_linear_algebra!(matrix::MacaulayMatrix, basis::Basis)
    resize!(matrix.coeffs, matrix.nlow)
    exact_sparse_rref!(matrix, basis)
end

#------------------------------------------------------------------------------

function interreduce_matrix_rows!(matrix::MacaulayMatrix, basis::Basis)
    resize!(matrix.lowrows, matrix.ncols)
    resize!(matrix.up2coef, matrix.ncols)
    resize!(matrix.low2coef, matrix.ncols)
    resize!(matrix.coeffs, matrix.ncols)

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{Int}}(undef, matrix.ncols)
    @inbounds for i in 1:matrix.nrows
        pivs[matrix.uprows[i][1]] = matrix.uprows[i]
        matrix.low2coef[matrix.uprows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.up2coef[i]])
    end

    densecfs = Vector{UInt64}(undef, matrix.ncols)
    k = matrix.nrows

    for i in 1:matrix.ncols
        l = matrix.ncols - i + 1
        if isassigned(pivs, l)
            densecfs .= UInt64(0)
            cfs = matrix.coeffs[matrix.low2coef[l]]
            reducexps = pivs[l]
            startcol = reducexps[1]
            @inbounds for j in 1:length(reducexps)
                densecfs[reducexps[j]] = cfs[j]
            end

            zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecfs, matrix, basis, pivs, startcol, startcol)

            matrix.lowrows[k] = newrow
            matrix.coeffs[matrix.low2coef[l]] = newcfs
            pivs[l] = matrix.lowrows[k]
            k -= 1
        end
    end

    # hmm
    # TODO
    matrix.npivots = matrix.nrows
end


#------------------------------------------------------------------------------

function convert_hashes_to_columns!(
            matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # col2hash = matrix.col2hash
    hdata    = symbol_ht.hashdata
    load     = symbol_ht.load

    # monoms from symbolic table represent one column in the matrix

    col2hash = Vector{Int}(undef, load - 1)
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

    # sort columns
    # TODO
    sort_columns_by_hash!(col2hash, symbol_ht)

    matrix.nleft = k
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    for k in 1:length(col2hash)
        hdata[col2hash[k]].idx = k
    end

    nterms = 0
    for k in 1:matrix.nup
        row = matrix.uprows[k]

        @inbounds for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end

        # TODO: not needed for now
        nterms += length(row)
    end

    for k in 1:matrix.nlow
        row = matrix.lowrows[k]

        @inbounds for j in 1:length(row)
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

function exact_sparse_linear_algebra_isgroebner!(matrix, basis)
    ncols  = matrix.ncols
    nlow   = matrix.nlow
    nright = matrix.nright
    nleft  = matrix.nleft

    # known pivots
    # no_copy
    pivs = Vector{Vector{Int}}(undef, ncols)
    @inbounds for i in 1:matrix.nup
        pivs[i] = matrix.uprows[i]
    end

    # CHANGED in order to prevent bug
    # when several rows in the matrix are equal
    l2c_tmp = Vector{Int}(undef, max(ncols, matrix.nlow))
    @inbounds for i in 1:nlow
        l2c_tmp[matrix.lowrows[i][1]] = matrix.low2coef[i]
    end
    # CHANGED
    # no_copy
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots
    # (not discovered yet)
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

    densecoeffs = zeros(UInt64, ncols)

    #=
    @warn "before reducing low"
    dump(matrix, maxdepth=5)
    @warn "lowrow2coef" rowidx2coef
    =#

    for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref  = basis.coeffs[rowidx2coef[i]]

        k = 0

        # we load coefficients into dense array
        # into rowexps indices
        # TODO: move this
        densecoeffs .= UInt64(0)
        @inbounds for j in 1:length(rowexps)
            densecoeffs[rowexps[j]] = cfsref[j]
        end

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        # zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        # @warn "reduced " zeroed newrow newcfs
        # if fully reduced
        zeroed && continue

        return false
    end
    return true
end

#------------------------------------------------------------------------------

function exact_sparse_rref_nf!(
        matrix::MacaulayMatrix,
        tobereduced::Basis,
        basis::Basis)

    ncols  = matrix.ncols
    nlow   = matrix.nlow
    nright = matrix.nright
    nleft  = matrix.nleft

    # known pivots
    # no_copy
    pivs = Vector{Vector{Int}}(undef, ncols)
    @inbounds for i in 1:matrix.nup
        pivs[i] = matrix.uprows[i]
    end

    # CHANGED in order to prevent bug
    # when several rows in the matrix are equal
    l2c_tmp = Vector{Int}(undef, max(ncols, matrix.nlow))
    @inbounds for i in 1:nlow
        l2c_tmp[matrix.lowrows[i][1]] = matrix.low2coef[i]
    end
    # CHANGED
    # no_copy
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots
    # (not discovered yet)
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows

    densecoeffs = zeros(UInt64, ncols)

    #=
    @warn "before reducing low"
    dump(matrix, maxdepth=5)
    @warn "lowrow2coef" rowidx2coef
    =#

    for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref  = tobereduced.coeffs[rowidx2coef[i]]

        k = 0

        # we load coefficients into dense array
        # into rowexps indices
        # TODO: move this
        densecoeffs .= UInt64(0)
        @inbounds for j in 1:length(rowexps)
            densecoeffs[rowexps[j]] = cfsref[j]
        end

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        # zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1)
        # @warn "reduced " zeroed newrow newcfs
        # if fully reduced
        zeroed && continue

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        matrix.lowrows[i]  = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.low2coef[i] = i
    end
    matrix.npivots = matrix.nrows = matrix.size = matrix.nlow
end

function exact_sparse_linear_algebra_nf!(
        matrix::MacaulayMatrix,
        tobereduced::Basis,
        basis::Basis)

    resize!(matrix.coeffs, matrix.nlow)
    exact_sparse_rref_nf!(matrix, tobereduced, basis)
end

#------------------------------------------------------------------------------

function convert_matrix_rows_to_basis_elements!(
            matrix::MacaulayMatrix, basis::Basis,
            ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # we mutate basis array directly by adding new elements

    check_enlarge_basis!(basis, matrix.npivots)
    rows = matrix.lowrows
    crs = basis.ndone

    for i in 1:matrix.npivots
        @inbounds colidx = rows[i][1]

        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.col2hash)
        # TODO : a constant

        # an interesing way to find coefficients
        basis.coeffs[crs + i] = matrix.coeffs[matrix.low2coef[colidx]]
        basis.gens[crs + i] = matrix.lowrows[i]
    end

    basis.ntotal += matrix.npivots
end


function convert_matrix_rows_to_basis_elements_use_symbol!(
            matrix::MacaulayMatrix, basis::Basis)

    check_enlarge_basis!(basis, matrix.npivots)

    crs = basis.ndone
    rows = matrix.lowrows

    for i in 1:matrix.npivots
        row = rows[i]
        colidx = row[1]

        @inbounds for j in 1:length(row)
            row[j] = matrix.col2hash[row[j]]
        end

        basis.coeffs[crs + i] = matrix.coeffs[matrix.low2coef[colidx]]
        basis.gens[crs + i] = row
    end
end

#------------------------------------------------------------------------------

function convert_nf_rows_to_basis_elements!(
        matrix::MacaulayMatrix, basis::Basis,
        ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    check_enlarge_basis!(basis, matrix.npivots)

    for i in 1:matrix.npivots
        basis.ndone += 1
        basis.nlead += 1
        basis.nonred[basis.nlead] = basis.ndone
        if isassigned(matrix.coeffs, i)
            row = matrix.lowrows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.col2hash)

            colidx = row[1]
            basis.coeffs[basis.ndone] = matrix.coeffs[i]
            basis.gens[basis.ndone] = row
        else
            empty!(basis.coeffs[basis.ndone])
            empty!(basis.gens[basis.ndone])
        end
    end
end
