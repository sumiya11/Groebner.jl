

mutable struct DoubleMacaulayMatrix{C}
    # pivot -> row
    leftpivs::Vector{Vector{Int}}
    leftcoeffs::Vector{Vector{C}}

    # idx -> row
    rightrows::Vector{Vector{Int}}
    rightcoeffs::Vector{Vector{C}}

    # pivot -> idx
    pivot2idx::Vector{Int}

    #
    lefthash2col::Dict{Int, Int}
    leftcol2hash::Vector{Int}

    righthash2col::Dict{Int, Int}
    rightcol2hash::Vector{Int}

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


function initialize_double_matrix(basis::Basis{C}) where {C<:Coeff}
    n = length(basis.gens)
    leftrows   = Vector{Vector{Int}}(undef, n)
    leftcoeffs = Vector{Vector{C}}(undef, n)

    pivot2idx = Vector{Int}(undef, n)

    m = length(basis.gens)
    rightrows  = Vector{Vector{Int}}(undef, m)
    rightcoeffs= Vector{Vector{C}}(undef, m)

    lsz = n
    lefthash2col = Dict{Int, Int}()
    leftcol2hash = Vector{Int}(undef, lsz)

    rsz = m
    righthash2col = Dict{Int, Int}()
    rightcol2hash = Vector{Int}(undef, rsz)

    DoubleMacaulayMatrix(leftrows, leftcoeffs, rightrows,
                rightcoeffs, pivot2idx, lefthash2col, leftcol2hash,
                righthash2col, rightcol2hash, n, 0, m, 0, 0)
end


function convert_to_double_dense_row(matrix, monom, vector::Basis{C}, ht) where {C<:Coeff}

    if matrix.nrcols >= length(matrix.rightcol2hash)
        resize!(matrix.rightcol2hash, 2*length(matrix.rightcol2hash))
    end
    matrix.nrcols += 1
    matrix.rightcol2hash[matrix.nrcols] = monom
    matrix.righthash2col[monom] = matrix.nrcols

    rightrow = zeros(C, matrix.nrcols)
    rightrow[end] = one(C)

    exps, coeffs = vector.gens[1], vector.coeffs[1]
    for i in 1:length(exps)
        if !haskey(matrix.lefthash2col, exps[i])
            if matrix.nlcols >= length(matrix.leftcol2hash)
                resize!(matrix.leftcol2hash, 2*length(matrix.leftcol2hash))
            end
            matrix.nlcols += 1
            matrix.lefthash2col[exps[i]] = matrix.nlcols
            matrix.leftcol2hash[matrix.nlcols] = exps[i]
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
function reduce_by_pivot_simultaneous!(leftrow, leftexps, leftcfs::Vector{CoeffFF},
                rightrow, rightexps, rightcfs, magic)

    # mul = -densecoeffs[i]
    # actually.. not bad!
    mul = (magic.divisor - leftrow[leftexps[1]]) % magic

    @inbounds for j in 1:length(leftexps)
        idx = leftexps[j]
        leftrow[idx] = (leftrow[idx] + mul*leftcfs[j]) % magic
    end

    @inbounds for j in 1:length(rightexps)
        idx = rightexps[j]
        rightrow[idx] = (rightrow[idx] + mul*rightcfs[j]) % magic
    end

    mul
end

#
# Finite field magic specialization
function normalize_double_row_sparse!(leftcfs::Vector{CoeffFF}, rightcfs, magic)
    pinv = invmod(leftcfs[1], magic.divisor) % magic
    @inbounds for i in 2:length(leftcfs)
        # row[i] *= pinv
        leftcfs[i] = (leftcfs[i] * pinv) % magic
    end
    @inbounds leftcfs[1] = one(leftcfs[1])

    @inbounds for i in 1:length(rightcfs)
        # row[i] *= pinv
        rightcfs[i] = (rightcfs[i] * pinv) % magic
    end
end

#
# Finite field magic specialization
function normalize_double_row_sparse!(leftcfs::Vector{CoeffQQ}, rightcfs, magic)
    pinv = inv(leftcfs[1])
    @inbounds for i in 2:length(leftcfs)
        # row[i] *= pinv
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
function reduce_by_pivot_simultaneous!(leftrow, leftexps, leftcfs::Vector{CoeffQQ},
                rightrow, rightexps, rightcfs, magic)

    # mul = -densecoeffs[i]
    # actually.. not bad!

    mul = -leftrow[leftexps[1]]

    @inbounds for j in 1:length(leftexps)
        idx = leftexps[j]
        leftrow[idx] = leftrow[idx] + mul*leftcfs[j]
    end

    @inbounds for j in 1:length(rightexps)
        idx = rightexps[j]
        rightrow[idx] = rightrow[idx] + mul*rightcfs[j]
    end

    mul
end


function reduce_double_dense_row_by_known_pivots_sparse!(
            matrix::DoubleMacaulayMatrix{C},
            leftrow, rightrow, magic) where {C}

    leftrows  = matrix.leftpivs
    rightrows = matrix.rightrows

    pivot2idx = matrix.pivot2idx

    # new row nonzero elements count
    k = 0
    uzero = C(0)

    # new pivot index
    np = -1

    if debug()
        @warn "in reduce" matrix.nlcols matrix.nrcols matrix.leftpivs
        @warn "hmm" leftrow
    end

    for i in 1:matrix.nlcols

        # if row element zero - no reduction
        @inbounds if leftrow[i] == uzero
            continue
        end

        # TODO: check this first?
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

        mul = reduce_by_pivot_simultaneous!(leftrow, leftexps, leftcfs,
                        rightrow, rightexps, rightcfs, magic)

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
    newrow = Vector{Int}(undef, k)
    newcfs = Vector{C}(undef, k)

    # store new row in sparse format
    # where k - number of structural nonzeros in new reduced row, k > 0
    j = 1
    @inbounds for i in np:length(row) # from new pivot
        @inbounds if row[i] != 0
            newrow[j] = i
            newcfs[j] = row[i]
            j += 1
        end
    end

    newrow, newcfs, j - 1
end

function linear_relation!(
            matrix::DoubleMacaulayMatrix,
            monom::Int, vector::Basis{C},
            ht) where {C<:Coeff}

    magic = select_divisor(vector.coeffs, vector.ch)

    leftrow, rightrow = convert_to_double_dense_row(matrix, monom, vector, ht)

    if debug()
        @warn "start"
        println(monom)
        println(vector.gens, " ", vector.coeffs)
        println(leftrow)
        println(rightrow)
        println(matrix)
    end

    reduced, np, k = reduce_double_dense_row_by_known_pivots_sparse!(matrix, leftrow, rightrow, magic)

    if debug()
        @warn "reduced"
        println(reduced, " ", np, " ", k)
        println(leftrow)
        println(rightrow)
    end

    if reduced
        # pass
    else
        lexps, lcoeffs, _ = extract_sparse_row(leftrow, np, k)
        rexps, rcoeffs, _ = extract_sparse_row(rightrow)

        normalize_double_row_sparse!(lcoeffs, rcoeffs, magic)

        if debug()
            @warn "extracted"
            println(lexps, " ", lcoeffs)
            println(rexps, " ", rcoeffs)
        end

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

    return reduced, rightrow
end
