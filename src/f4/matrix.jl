
const ColumnIdx = Int32

mutable struct MacaulayMatrix{T <: Coeff}
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

    up2mult::Vector{MonomIdx}
    low2mult::Vector{MonomIdx}
end

# Initializes an empty matrix with coefficients of type T
function initialize_matrix(ring::PolyRing, ::Type{T}) where {T <: Coeff}
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

    up2mult = Vector{MonomIdx}(undef, 0)
    low2mult = Vector{MonomIdx}(undef, 0)

    MacaulayMatrix(
        uprows,
        lowrows,
        col2hash,
        coeffs,
        size,
        npivots,
        nrows,
        ncols,
        nup,
        nlow,
        nleft,
        nright,
        up2coef,
        low2coef,
        up2mult,
        low2mult
    )
end

# Refresh and initialize matrix for `npairs` elements
function reinitialize_matrix!(matrix::MacaulayMatrix{T}, npairs::Int) where {T}
    resize!(matrix.uprows, npairs * 2)
    resize!(matrix.lowrows, npairs * 2)
    resize!(matrix.up2coef, npairs * 2)
    resize!(matrix.low2coef, npairs * 2)
    resize!(matrix.up2mult, npairs * 2)
    resize!(matrix.low2mult, npairs * 2)
    matrix.size = 2 * npairs
    matrix.ncols = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nup = 0
    matrix.nlow = 0
    matrix
end

function linear_algebra!(
    graph,
    ring::PolyRing,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg,
    rng
)
    resize!(matrix.coeffs, matrix.nlow)
    if linalg === :learn
        return learn_sparse_rref!(graph, ring, matrix, basis, rng)
    else
        return apply_sparse_rref!(graph, ring, matrix, basis, rng)
    end
    true
end

function linear_algebra!(ring::PolyRing, matrix::MacaulayMatrix, basis::Basis, linalg, rng)
    resize!(matrix.coeffs, matrix.nlow)
    if linalg === :deterministic
        exact_sparse_rref!(ring, matrix, basis)
        return true
    else
        randomized_sparse_rref!(ring, matrix, basis, rng)
        return true
    end
    true
end

#------------------------------------------------------------------------------

# Normalize `row` by first coefficient
#
# Finite field magic specialization
function normalize_sparse_row!(row::Vector{T}, arithmetic) where {T <: CoeffFF}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = mod_x(invmod(row[1], divisor(arithmetic)), arithmetic)
    @inbounds for i in 2:length(row)
        row[i] = mod_x(row[i] * pinv, arithmetic)
    end
    @inbounds row[1] = one(row[1])
    row
end

# Normalize `row` by first coefficient
#
# Rational field specialization
function normalize_sparse_row!(row::Vector{T}, ch) where {T <: CoeffQQ}
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

# reduces row by mul*cfs modulo ch at indices positions
#
# Finite field magic specialization
function reduce_by_pivot!(
    row::Vector{T},
    indices::Vector{ColumnIdx},
    cfs::Vector{T},
    arithmetic
) where {T <: CoeffFF}
    @inbounds mul = divisor(arithmetic) - row[indices[1]]

    # on our benchmarks usually
    # length(row) / length(indices) varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = mod_x(row[idx] + mul * cfs[j], arithmetic)
    end

    nothing
end

# reduces row by mul*cfs at indices positions
#
# Rational field specialization
function reduce_by_pivot!(
    row::Vector{T},
    indices::Vector{ColumnIdx},
    cfs::Vector{T},
    arithmetic::BuiltinArithmeticQQ
) where {T <: CoeffQQ}
    @inbounds mul = -row[indices[1]]
    buf2 = arithmetic.buf2
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * cfs[j]
        # Base.GMP.MPQ.mul!(buf2, cfs[j], mul)
        # row[idx] = row[idx] + buf2
    end
    nothing
end

# zero entries of densecoeffs and load coefficients cfsref to indices rowexps
#
# Finite field specialization.
function load_indexed_coefficients!(
    densecoeffs::Vector{T},
    rowexps,
    cfsref
) where {T <: CoeffFF}
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
function load_indexed_coefficients!(
    densecoeffs::Vector{T},
    rowexps,
    cfsref
) where {T <: CoeffQQ}
    densecoeffs .= T(0)
    @inbounds for j in 1:length(rowexps)
        densecoeffs[rowexps[j]] = cfsref[j]
    end
end

# 
function extract_sparse_row!(inds, vals, dense::Vector{T}, from, to) where {T}
    z = zero(T)
    j = 1
    @inbounds for i in from:to # starting from new pivot
        if dense[i] != z
            inds[j] = i
            vals[j] = dense[i]
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

function reduce_dense_row_by_known_pivots_sparse!(
    densecoeffs::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs::Vector{Vector{ColumnIdx}},
    startcol,
    tmp_pos,
    magic;
    exact_colmap::Bool=false
) where {C <: Coeff}
    reduce_dense_row_by_known_pivots_sparse!(
        densecoeffs,
        matrix,
        basis,
        pivs,
        ColumnIdx(startcol),
        ColumnIdx(tmp_pos),
        magic;
        exact_colmap=exact_colmap
    )
end

function reduce_dense_row_by_known_pivots_sparse!(
    reducers,
    densecoeffs::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs::Vector{Vector{ColumnIdx}},
    startcol,
    tmp_pos,
    arithmetic;
    exact_colmap::Bool=false
) where {C <: Coeff}
    ncols = matrix.ncols
    nleft = matrix.nleft

    # new row nonzero elements count
    k = 0
    uzero = zero(C)

    # new pivot index
    np = -1

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
        elseif i <= nleft # if reducer is from the upper part of the matrix
            @inbounds cfs = basis.coeffs[matrix.up2coef[i]]
            push!(reducers, (matrix.up2coef[i], matrix.up2mult[i]))
        else # if reducer is from the lower part of the matrix
            @inbounds cfs = matrix.coeffs[matrix.low2coef[i]]
        end

        reduce_by_pivot!(densecoeffs, reducerexps, cfs, arithmetic)
    end

    newrow = Vector{ColumnIdx}(undef, k)
    newcfs = Vector{C}(undef, k)

    # all reduced !
    if k == 0
        return true, newrow, newcfs
    end

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    @assert k > 0
    extract_sparse_row!(newrow, newcfs, densecoeffs, startcol, ncols)
    return false, newrow, newcfs
end

function reduce_dense_row_by_known_pivots_sparse!(
    densecoeffs::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs::Vector{Vector{ColumnIdx}},
    startcol::ColumnIdx,
    tmp_pos::ColumnIdx,
    arithmetic;
    exact_colmap::Bool=false
) where {C <: Coeff}
    ncols = matrix.ncols
    nleft = matrix.nleft

    # new row nonzero elements count
    k = 0
    uzero = zero(C)

    # new pivot index
    np = -1

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
        elseif i <= nleft # if reducer is from the upper part of the matrix
            @inbounds cfs = basis.coeffs[matrix.up2coef[i]]
        else # if reducer is from the lower part of the matrix
            @inbounds cfs = matrix.coeffs[matrix.low2coef[i]]
        end

        reduce_by_pivot!(densecoeffs, reducerexps, cfs, arithmetic)
    end

    newrow = Vector{ColumnIdx}(undef, k)
    newcfs = Vector{C}(undef, k)

    # all reduced !
    if k == 0
        return true, newrow, newcfs
    end

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    @assert k > 0
    extract_sparse_row!(newrow, newcfs, densecoeffs, startcol, ncols)
    return false, newrow, newcfs
end

#------------------------------------------------------------------------------

function absolute_pivots!(matrix::MacaulayMatrix)
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:(matrix.nup)
        pivs[i] = matrix.uprows[i]
    end
    l2c_tmp = Vector{ColumnIdx}(undef, max(matrix.ncols, matrix.nlow))
    @inbounds for i in 1:(matrix.nlow)
        l2c_tmp[matrix.lowrows[i][1]] = matrix.low2coef[i]
    end
    pivs, l2c_tmp
end

function interreduce_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs,
    arithmetic;
    reversed=false
) where {C <: Coeff}
    # number of new pivots
    newpivs = 0
    densecoeffs = zeros(C, matrix.ncols)
    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)
    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:(matrix.nright)
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

        _, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            startcol,
            arithmetic
        )
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

function interreduce_lower_part_apply!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs,
    arithmetic;
    reversed=false
) where {C <: Coeff}
    # number of new pivots
    newpivs = 0
    densecoeffs = zeros(C, matrix.ncols)
    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)
    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:(matrix.nright)
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

        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            startcol,
            arithmetic
        )
        newpivs += 1
        if zeroed
            # unlucky cancellation
            return false
        end

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
    true
end

function interreduce_lower_part_learn!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivs,
    arithmetic;
    reversed=false
) where {C <: Coeff}
    # number of new pivots
    newpivs = 0
    densecoeffs = zeros(C, matrix.ncols)
    # a row to be reduced for each column
    resize!(matrix.lowrows, matrix.nright)

    not_reduced_to_zero = Int[]

    # interreduce new pivots..
    # .. for each right (non-pivotal) column
    @inbounds for i in 1:(matrix.nright)
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

        _, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            startcol,
            arithmetic
        )
        newpivs += 1

        push!(not_reduced_to_zero, k)

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

    push!(graph.matrix_infos, (nup=matrix.nup, nlow=matrix.nlow, ncols=matrix.ncols))
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (matrix.up2coef[1:(matrix.nup)], matrix.up2mult[1:(matrix.nup)])
    )
    push!(graph.matrix_lower_rows, (Int[], Int[]))

    # shrink matrix
    matrix.npivots = matrix.nrows = matrix.size = newpivs
    resize!(matrix.lowrows, newpivs)
end

# Linear algebra option 1
function exact_sparse_rref!(
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    ncols = matrix.ncols
    nlow = matrix.nlow

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

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
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
            arithmetic
        )
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
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
end

# Linear algebra option 3
function learn_sparse_rref!(
    graph::ComputationGraphF4,
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    rng
) where {C <: Coeff}
    ncols = matrix.ncols
    nlow = matrix.nlow

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = absolute_pivots!(matrix)
    @log level = -6 "absolute_pivots!" pivs l2c_tmp

    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lowrows
    densecoeffs = zeros(C, ncols)

    not_reduced_to_zero = Int[]
    useful_reducers = Set{Int}()

    @log level = -6 "Low to coef" rowidx2coef matrix.low2mult

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
        load_indexed_coefficients!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        reducers = Any[]
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            reducers,
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
            arithmetic
        )
        # if fully reduced
        zeroed && continue

        # NOTE: we are not adding reducers from lowrows!

        push!(not_reduced_to_zero, i)
        for rr in reducers
            push!(useful_reducers, rr)
        end

        @log level = -7 "Not zero" i rowidx2coef[i] matrix.low2mult[i]

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.low2coef[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    useful_reducers_sorted = sort(collect(useful_reducers))
    @log level = -7 "" useful_reducers_sorted
    push!(graph.matrix_infos, (nup=matrix.nup, nlow=matrix.nlow, ncols=matrix.ncols))
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (map(first, useful_reducers_sorted), map(last, useful_reducers_sorted))
    )
    push!(
        graph.matrix_lower_rows,
        (rowidx2coef[not_reduced_to_zero], matrix.low2mult[not_reduced_to_zero])
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
    ncols = matrix.ncols
    nlow = matrix.nlow

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    # YES
    pivs, l2c_tmp = absolute_pivots!(matrix)
    @log level = -7 "absolute_pivots!" pivs l2c_tmp

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

        @log level = -7 "$i from $nlow" rowexps cfsref rowidx2coef

        # we load coefficients into dense array
        # into rowexps indices
        load_indexed_coefficients!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.uprows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
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
        matrix.low2coef[newrow[1]] = i

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
    ncols = matrix.ncols
    nlow = matrix.nlow

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

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

            # reduce it with known pivots from matrix.uprows
            # first nonzero in densecoeffs is at startcol position
            zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
                densecoeffs,
                matrix,
                basis,
                pivs,
                startcol,
                -1,
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
            matrix.low2coef[newrow[1]] = absolute_i

            # normalize if needed
            normalize_sparse_row!(matrix.coeffs[absolute_i], arithmetic)

            ctr += 1
        end
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
end

function exact_sparse_rref_interreduce_apply!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lowrows, matrix.ncols)
    resize!(matrix.up2coef, matrix.ncols)
    resize!(matrix.up2mult, matrix.ncols)
    resize!(matrix.low2coef, matrix.ncols)
    resize!(matrix.low2mult, matrix.ncols)
    resize!(matrix.coeffs, matrix.ncols)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.uprows[i][1]] = matrix.uprows[i]
        matrix.low2coef[matrix.uprows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.up2coef[i]])
    end

    interreduce_lower_part_apply!(matrix, basis, pivs, arithmetic, reversed=true)
end

function exact_sparse_rref_interreduce!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lowrows, matrix.ncols)
    resize!(matrix.up2coef, matrix.ncols)
    resize!(matrix.up2mult, matrix.ncols)
    resize!(matrix.low2coef, matrix.ncols)
    resize!(matrix.low2mult, matrix.ncols)
    resize!(matrix.coeffs, matrix.ncols)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.uprows[i][1]] = matrix.uprows[i]
        matrix.low2coef[matrix.uprows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.up2coef[i]])
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic, reversed=true)
end

function exact_sparse_rref_interreduce_learn!(
    graph,
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.lowrows, matrix.ncols)
    resize!(matrix.up2coef, matrix.ncols)
    resize!(matrix.up2mult, matrix.ncols)
    resize!(matrix.low2coef, matrix.ncols)
    resize!(matrix.low2mult, matrix.ncols)
    resize!(matrix.coeffs, matrix.ncols)

    arithmetic = select_arithmetic(matrix.coeffs, UInt64(ring.ch))

    # same pivs as for rref
    # pivs: column idx --> vector of present columns
    pivs = Vector{Vector{ColumnIdx}}(undef, matrix.ncols)
    @inbounds for i in 1:(matrix.nrows)
        pivs[matrix.uprows[i][1]] = matrix.uprows[i]
        matrix.low2coef[matrix.uprows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.up2coef[i]])
    end

    interreduce_lower_part_learn!(graph, matrix, basis, pivs, arithmetic, reversed=true)
end

function exact_sparse_rref_isgroebner!(
    ring::PolyRing,
    matrix::MacaulayMatrix{C},
    basis::Basis{C}
) where {C <: Coeff}
    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    pivs, l2c_tmp = absolute_pivots!(matrix)
    rowidx2coef = matrix.low2coef
    matrix.low2coef = l2c_tmp

    upivs = matrix.lowrows

    densecoeffs = zeros(C, matrix.ncols)

    @inbounds for i in 1:(matrix.nlow)
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
        zeroed, _, _ = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
            arithmetic
        )
        # # @warn "reduced " zeroed newrow newcfs
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
    basis::Basis{C}
) where {C <: Coeff}
    resize!(matrix.coeffs, matrix.nlow)

    ncols = matrix.ncols
    nlow = matrix.nlow

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

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
        zeroed, newrow, newcfs = reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
            arithmetic
        )
        # # @warn "reduced " zeroed
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

function column_to_monom_mapping!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    col2hash = Vector{MonomIdx}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        # column to hash index
        col2hash[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == 2
            k += 1
        end
    end

    sort_columns_by_hash!(col2hash, symbol_ht)

    matrix.nleft = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(col2hash)
        hdata[col2hash[k]].idx = k
    end

    @inbounds for k in 1:(matrix.nup)
        row = matrix.uprows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nlow)
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
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)

    # we mutate basis array directly by adding new elements

    resize_basis_if_needed!(basis, matrix.npivots)
    rows = matrix.lowrows
    crs = basis.nprocessed

    @inbounds for i in 1:(matrix.npivots)
        colidx = rows[i][1]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.col2hash)
        basis.coeffs[crs + i] = matrix.coeffs[matrix.low2coef[colidx]]
        basis.monoms[crs + i] = matrix.lowrows[i]
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
            row = matrix.lowrows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.col2hash)
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
    row::Vector{ColumnIdx},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    col2hash::Vector{MonomIdx}
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
        hidx = col2hash[row[l]]

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
