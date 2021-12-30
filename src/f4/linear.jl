
#------------------------------------------------------------------------------

mutable struct Matrix
    #=
        Matrix of the following structure

        | A  B |
        | C  D |

        where A consists of pivots corresponding to
        first leading monomials in reduction
    =#

    size::Int
    npivots::Int
    nrows::Int
    ncols::Int

    # number of upper rows (in A and B section)
    nup::Int
    # number of lower rows (in C and D section)
    nlow::Int
    # number of left cols  (in A and C section)
    nleft::Int
    # number of right cols (in B and D section)
    nright::Int

    reducer
    reduced

    col2hash

    coeffs # hmm

end

function initialize_matrix()
    return Matrix(0, 0, 0, 0, 0, 0, 0, 0, [], [], [], [])
end

function reinitialize_matrix!(matrix, npairs)
    resize!(matrix.reducer, npairs*2)
    resize!(matrix.reduced, npairs*2)
    matrix.size = 2*npairs
    matrix
end

#------------------------------------------------------------------------------

function normalize_sparse_row!(row)
    pinv = inv(row[1])
    for i in 2:length(row)
        row[i] *= pinv
    end
    row[1] = 1
    row
end

function reduce_dense_row_by_known_pivots_sparse!(
            dr, matrix, basis, pivs, dpiv, tmp_pos)

    ncols = matrix.ncols
    nleft = matrix.nleft
    k = 1
    np = -1

    for i in dpiv:ncols
        # if zero - no reduction
        if dr[i] == 0
            continue
        end
        # if pivot not defined
        if pivs[i] == 0
            if np == -1
                np = i
            end
            k += 1
            continue
        end

        mul = inv(dr[i])
        dts = pivs[i]

        if i <= nleft
            cfs = basis.coeffs[dts[1]]
        else
            cfs = matrix.coeffs
        end

        for j in 1:length(dts)
            dr[ds[j]] += mul * cfs[j]
        end
    end

    if k == 0
        return C_NULL
    end

    # store new row in sparse format
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
    pivs = copy(matrix.reducer)

    # unknown pivots
    # (not discovered yet)
    upivs = matrix.reduced

    col = 1

    for i in 1:nlow
        @debug "low row $i.."

        npiv = upivs[i]
        # corresponding coefficients from basis
        # we need to locate it somehow
        cfs  = basis.coeffs[i] # TODO

        k = 0
        drl = zeros(ground, ncols)

        for j in 1:length(npiv)
            drl[npiv[j]] = cfs[j]
        end

        npiv = reduce_dense_row_by_known_pivots_sparse!(drl, matrix, basis, pivs, col, i)
        (npiv == 0) && break

        if matrix.coeffs[npiv[1]][1] != 1
            normalize_sparse_row!(matrix.coeffs[npiv[1]])
        end

        cfs = matrix.coeffs[npiv[1]]
    end

    # number of new pivots
    npivs = 0

    # interreduce new pivots
    for i in 1:nright
        k = ncols - 1 - i
        if pivs[k]
            col = ds[1]
            dr = zeros(Int, ncols)
            cfs = matrix.coeffs[pivs[k][1]]
            for j in 1:length(cfs)
                dr[ds[j]] = cfs[j]
            end
            npivs += 1
            matrix.reduced[npivs] = reduce_dense_row_by_known_pivots_sparse!(drl, matrix, basis, pivs, col, 1)
        end
    end

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = npivs
    resize!(matrix.reduced, npivs)
end

function exact_sparse_linear_algebra!(matrix, basis)
    resize!(matrix.reduced, matrix.nlow)
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
        row = matrix.reducer[k]

        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end

        # TODO
        nterms += length(row)
    end

    for k in 1:matrix.nlow
        row = matrix.reduced[k]

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


function convert_matrix_rows_to_basis_elements(
            matrix, state, basis_ht, symbol_ht, col2hash)

    # we mutate basis array directly by adding new elements

    crs = state.basis.load

    for i in 1:matrix.np
        insert_in_basis_hash_table_pivots(rows[i], basis_ht, symbol_ht, col2hash)
        # TODO : a constant

        state.basis.coeffs[crs + i]      = matrix.coeffs[i]
        state.basis.basis_polys[crs + i] = matrix.reduced[i]
    end
end
