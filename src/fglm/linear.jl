

mutable struct DoubleMacaulayMatrix
    leftrows::Vector{Vector{Int}}
    rightrows::Vector{Vector{Int}}

    leftcoeffs::Vector{Vector{UInt64}}
    rightcoeffs::Vector{Vector{UInt64}}

    lefthash2col::Vector{Int}
    leftcol2hash::Vector{Int}

    righthash2col::Vector{Int}
    rightcol2hash::Vector{Int}

    leftsize::Int
    rightsize::Int
end


function initialize_double_matrix(ring::PolyRing)
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

    DoubleMacaulayMatrix(uprows, lowrows, col2hash, coeffs,
            size, npivots, nrows, ncols,
            nup, nlow, nleft, nright,
            up2coef, low2coef)
end


function linear_relation!(
            matrix::DoubleMacaulayMatrix,
            monom::Int, vector::Basis)
    denserow = convert_vector_to_dense_row(matrix, monom, vector)
    reduce_double_dense_row_by_known_pivots_sparse!(matrix, denserow)

    if there is a relation
        return true
    else
        false
end
