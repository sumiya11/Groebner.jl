# This file is a part of Groebner.jl. License is GNU GPL v2.

# Left part -- linear combinations of residuals (normal forms)
# Right part -- linear combinations of monomials (preimages)
mutable struct WideMacaulayMatrix{C}
    leftpivs::Vector{Vector{Int}}
    leftcoeffs::Vector{Vector{C}}

    rightrows::Vector{Vector{Int}}
    rightcoeffs::Vector{Vector{C}}

    pivot2idx::Vector{Int}

    lefthash2col::Dict{Int, Int}
    leftcolumn_to_monom::Vector{Int}

    righthash2col::Dict{Int, Int}
    rightcolumn_to_monom::Vector{Int}

    # size of leftpivs
    nlsize::Int

    # load of rightrows
    nrrows::Int
    # size of rightrows
    nrsize::Int

    # left cols filled
    nlcols::Int
    # right cols filled
    nrcols::Int
end

function wide_matrix_initialize(basis::Basis{C}) where {C <: Coeff}
    n          = length(basis.monoms)
    leftrows   = Vector{Vector{Int}}(undef, n)
    leftcoeffs = Vector{Vector{C}}(undef, n)

    pivot2idx = Vector{Int}(undef, n)

    m           = length(basis.monoms)
    rightrows   = Vector{Vector{Int}}(undef, m)
    rightcoeffs = Vector{Vector{C}}(undef, m)

    lsz = n
    lefthash2col = Dict{Int, Int}()
    leftcolumn_to_monom = Vector{Int}(undef, lsz)

    rsz = m
    righthash2col = Dict{Int, Int}()
    rightcolumn_to_monom = Vector{Int}(undef, rsz)

    WideMacaulayMatrix(
        leftrows,
        leftcoeffs,
        rightrows,
        rightcoeffs,
        pivot2idx,
        lefthash2col,
        leftcolumn_to_monom,
        righthash2col,
        rightcolumn_to_monom,
        n,
        0,
        m,
        0,
        0
    )
end

function convert_to_wide_dense_row(matrix, monom, vector::Basis{C}) where {C <: Coeff}
    if matrix.nrcols >= length(matrix.rightcolumn_to_monom)
        resize!(matrix.rightcolumn_to_monom, 2 * length(matrix.rightcolumn_to_monom))
    end
    matrix.nrcols += 1
    matrix.rightcolumn_to_monom[matrix.nrcols] = monom
    matrix.righthash2col[monom] = matrix.nrcols

    rightrow = zeros(C, matrix.nrcols)
    rightrow[end] = one(C)

    exps, coeffs = vector.monoms[1], vector.coeffs[1]
    for i in 1:length(exps)
        if !haskey(matrix.lefthash2col, exps[i])
            if matrix.nlcols >= length(matrix.leftcolumn_to_monom)
                resize!(matrix.leftcolumn_to_monom, 2 * length(matrix.leftcolumn_to_monom))
            end
            matrix.nlcols += 1
            matrix.lefthash2col[exps[i]] = matrix.nlcols
            matrix.leftcolumn_to_monom[matrix.nlcols] = exps[i]
        end
    end

    leftrow = zeros(C, matrix.nlcols)
    for i in 1:length(exps)
        leftrow[matrix.lefthash2col[exps[i]]] = coeffs[i]
    end

    leftrow, rightrow
end

# reduces row by mul*cfs modulo ch at indices positions
#
# Finite field magic specialization
function reduce_by_pivot_simultaneous!(
    leftrow,
    leftexps,
    leftcfs::Vector{T},
    rightrow,
    rightexps,
    rightcfs,
    arithmetic
) where {T <: CoeffZp}
    mul = mod_p(divisor(arithmetic) - leftrow[leftexps[1]], arithmetic)

    @inbounds for j in 1:length(leftexps)
        idx = leftexps[j]
        leftrow[idx] = mod_p(leftrow[idx] + mul * leftcfs[j], arithmetic)
    end

    @inbounds for j in 1:length(rightexps)
        idx = rightexps[j]
        rightrow[idx] = mod_p(rightrow[idx] + mul * rightcfs[j], arithmetic)
    end

    mul
end

function normalize_wide_row_sparse!(
    leftcfs::Vector{T},
    rightcfs,
    arithmetic
) where {T <: CoeffZp}
    pinv = mod_p(invmod(leftcfs[1], divisor(arithmetic)), arithmetic)
    @inbounds for i in 2:length(leftcfs)
        leftcfs[i] = mod_p(leftcfs[i] * pinv, arithmetic)
    end
    @inbounds leftcfs[1] = one(leftcfs[1])

    @inbounds for i in 1:length(rightcfs)
        rightcfs[i] = mod_p(rightcfs[i] * pinv, arithmetic)
    end
end

function normalize_wide_row_sparse!(
    leftcfs::Vector{T},
    rightcfs,
    arithmetic
) where {T <: CoeffQQ}
    pinv = inv(leftcfs[1])
    @inbounds for i in 2:length(leftcfs)
        leftcfs[i] = leftcfs[i] * pinv
    end
    @inbounds leftcfs[1] = one(leftcfs[1])

    @inbounds for i in 1:length(rightcfs)
        rightcfs[i] = rightcfs[i] * pinv
    end
end

# reduces row by mul*cfs modulo ch at indices positions
#
# Rational field specialization
function reduce_by_pivot_simultaneous!(
    leftrow,
    leftexps,
    leftcfs::Vector{T},
    rightrow,
    rightexps,
    rightcfs,
    arithmetic
) where {T <: CoeffQQ}
    mul = -leftrow[leftexps[1]]

    @inbounds for j in 1:length(leftexps)
        idx = leftexps[j]
        leftrow[idx] = leftrow[idx] + mul * leftcfs[j]
    end

    @inbounds for j in 1:length(rightexps)
        idx = rightexps[j]
        rightrow[idx] = rightrow[idx] + mul * rightcfs[j]
    end

    mul
end

function reduce_wide_dense_row_by_known_pivots_sparse!(
    matrix::WideMacaulayMatrix{C},
    leftrow,
    rightrow,
    arithmetic
) where {C}
    leftrows  = matrix.leftpivs
    rightrows = matrix.rightrows

    pivot2idx = matrix.pivot2idx

    # new row nonzero elements count
    k = 0
    uzero = C(0)

    # new pivot index
    np = -1

    for i in 1:(matrix.nlcols)

        # if row element zero - no reduction
        @inbounds if leftrow[i] == uzero
            continue
        end

        # TODO: check this before checking leftrow[i] == uzero?
        if !isassigned(leftrows, i)
            if np == -1
                np = i
            end
            k += 1
            continue
        end

        # exponents of reducer row at column i
        leftexps = leftrows[i]
        leftcfs  = matrix.leftcoeffs[i]

        # here map pivot --> when added
        rightexps = rightrows[pivot2idx[i]]
        rightcfs  = matrix.rightcoeffs[pivot2idx[i]]

        reduce_by_pivot_simultaneous!(
            leftrow,
            leftexps,
            leftcfs,
            rightrow,
            rightexps,
            rightcfs,
            arithmetic
        )
    end

    return k == 0, np, k
end

function extract_sparse_row(row)
    newrow, newcfs, k = extract_sparse_row(row, 1, length(row))
    resize!(newrow, k)
    resize!(newcfs, k)
    newrow, newcfs, k
end

function extract_sparse_row(row::Vector{C}, np, k) where {C}
    newrow = Vector{MonomId}(undef, k)
    newcfs = Vector{C}(undef, k)

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    j = 1
    @inbounds for i in np:length(row) # from new pivot
        if row[i] != 0
            newrow[j] = i
            newcfs[j] = row[i]
            j += 1
        end
    end

    newrow, newcfs, j - 1
end

const _sparisty_factors = []

@timeit function find_linear_relation!(
    matrix::WideMacaulayMatrix,
    monom::MonomId,
    vector::Basis{C},
    arithmetic
) where {C <: Coeff}
    leftrow, rightrow = convert_to_wide_dense_row(matrix, monom, vector)

    reduced, np, k =
        reduce_wide_dense_row_by_known_pivots_sparse!(matrix, leftrow, rightrow, arithmetic)

    if reduced
        # pass
    else
        lexps, lcoeffs, _ = extract_sparse_row(leftrow, np, k)
        rexps, rcoeffs, _ = extract_sparse_row(rightrow)

        normalize_wide_row_sparse!(lcoeffs, rcoeffs, arithmetic)

        while np >= matrix.nlsize
            matrix.nlsize *= 2
            resize!(matrix.leftpivs, matrix.nlsize)
            resize!(matrix.leftcoeffs, matrix.nlsize)
            resize!(matrix.pivot2idx, matrix.nlsize)
        end
        matrix.leftpivs[np] = lexps
        matrix.leftcoeffs[np] = lcoeffs

        matrix.nrrows += 1
        matrix.pivot2idx[np] = matrix.nrrows

        if matrix.nrrows >= matrix.nrsize
            matrix.nrsize *= 2
            resize!(matrix.rightrows, matrix.nrsize)
            resize!(matrix.rightcoeffs, matrix.nrsize)
        end
        matrix.rightrows[matrix.nrrows] = rexps
        matrix.rightcoeffs[matrix.nrrows] = rcoeffs
    end

    nnz_left = 0
    nnz_left_rows = 0
    for i in 1:(matrix.nlsize)
        if isassigned(matrix.leftpivs, i)
            nnz_left += length(matrix.leftpivs[i])
            nnz_left_rows += 1
        end
    end
    nnz_right = 0
    nnz_right_rows = 0
    for i in 1:(matrix.nrsize)
        if isassigned(matrix.rightrows, i)
            nnz_right += length(matrix.rightrows[i])
            nnz_right_rows += 1
        end
    end
    dimsleft = (nnz_left_rows, matrix.nlcols)
    dimsright = (nnz_right_rows, matrix.nrcols)
    spleft = nnz_left / prod(dimsleft)
    spright = nnz_right / prod(dimsright)

    push!(
        _sparisty_factors,
        (
            nnz_left=nnz_left,
            dimsleft=dimsleft,
            spleft=spleft,
            nnz_right=nnz_right,
            dimsright=dimsright,
            spright=spright
        )
    )

    return reduced, rightrow
end

function extract_linear_basis(ring, matrix::WideMacaulayMatrix{C}) where {C}
    exps = Vector{Vector{MonomId}}(undef, matrix.nrrows)
    coeffs = Vector{Vector{C}}(undef, matrix.nrrows)

    for i in 1:(matrix.nrrows)
        exps[i] = matrix.rightrows[i]
        coeffs[i] = matrix.rightcoeffs[i]
        for j in 1:length(exps[i])
            exps[i][j] = matrix.rightcolumn_to_monom[exps[i][j]]
        end
    end

    linbasis = basis_initialize(ring, exps, coeffs)

    linbasis.nprocessed = length(exps)
    linbasis.nnonredundant = length(exps)
    linbasis.nonredundant = collect(1:(linbasis.nprocessed))

    linbasis
end
