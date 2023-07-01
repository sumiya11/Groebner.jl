# Sorting monomials, polynomials, and other things.

# Unstable algorithm should be a bit faster for large arrays
# (basically, quicksort vs. mergesort)
_default_sorting_alg() = Base.Sort.DEFAULT_UNSTABLE

# Sorts vector v at indices from:to
function sort_part!(
    v,
    from::Integer,
    to::Integer;
    lt=isless,
    alg=_default_sorting_alg(),
    by=identity,
    rev=false
)
    ordr = Base.Sort.ord(lt, by, rev, Base.Sort.Forward)
    sort!(v, from, to, alg, ordr)
end

#------------------------------------------------------------------------------

# Sorts polynomials from the basis
# by their leading monomial in the non-decreasing way
# by the given term ordering.
# Also sorts any arrays passed in `abc` in the same order as basis.
function sort_polys_by_lead_increasing!(
    basis::Basis,
    ht::MonomialHashtable,
    abc...;
    ord::Ord=ht.ord
) where {Ord <: AbstractMonomialOrdering}
    gens = basis.monoms
    exps = ht.monoms

    inds = collect(1:(basis.ntotal))

    cmps =
        (x, y) ->
            monom_isless(@inbounds(exps[gens[x][1]]), @inbounds(exps[gens[y][1]]), ord)

    sort!(inds, lt=cmps, alg=_default_sorting_alg())

    # use array assignment insted of elemewise assignment
    basis.monoms[1:(basis.ntotal)] = basis.monoms[inds]
    basis.coeffs[1:(basis.ntotal)] = basis.coeffs[inds]
    for a in abc
        a[1:(basis.ntotal)] = a[inds]
    end
    nothing
end

#------------------------------------------------------------------------------

# Sorts pairs from pairset in the range [from, from+sz]
# by lcm total degrees in increasing order
#
# Used in update once per one f4 iteration to sort pairs in pairset;
# Also used in normal selection strategy also once per iteration
function sort_pairset_by_degree!(ps::Pairset, from::Int, sz::Int)
    sort_part!(ps.pairs, from, from + sz, by=p -> p.deg, alg=_default_sorting_alg())
end

#------------------------------------------------------------------------------

# Sorts the first `npairs` pairs from `pairset` by non-decreasing order of
# the exponent vector of the lcm wrt. the given monomial ordering
function sort_pairset_by_lcm!(pairset::Pairset, npairs::Int, ht::MonomialHashtable)
    exps = ht.exponents

    cmps = (x, y) -> monom_isless(@inbounds(exps[x.lcm]), @inbounds(exps[y.lcm]), ht.ord)

    sort_part!(pairset.pairs, 1, npairs, lt=cmps, alg=_default_sorting_alg())
end

#------------------------------------------------------------------------------

# Sorts generators selected in normal strategy function
# by their ordering in the current basis (identity sort)
function sort_generators_by_position!(gens::Vector{Int}, load::Int)
    sort_part!(gens, 1, load, alg=_default_sorting_alg())
end

#------------------------------------------------------------------------------
# Sorting matrix rows and columns.
# See f4/matrix.jl for details.

# Comparator for matrix rows a and b.
# A row matrix is encoded by the sparse array of 
# its nonzero column indices.
function matrix_row_decreasing_cmp(a::Vector{T}, b::Vector{T}) where {T <: ColumnIdx}
    #= a, b - rows as arrays of exponent indices =#
    # va and vb are the leading columns
    @inbounds va = a[1]
    @inbounds vb = b[1]

    if va > vb
        return false
    end
    if va < vb
        return true
    end

    # if same column index => compare the density of rows
    va = length(a)
    vb = length(b)
    if va > vb
        return true
    end
    if va < vb
        return false
    end

    # hmm, equal rows?
    # should never happen
    return false
end

# Comparator for matrix rows a and b.
# A row matrix is encoded by the sparse array of 
# its nonzero column indices.
function matrix_row_increasing_cmp(a::Vector{T}, b::Vector{T}) where {T <: ColumnIdx}
    #= a, b - rows as arrays of exponent hashes =#
    # va and vb are the leading columns
    @inbounds va = a[1]
    @inbounds vb = b[1]

    if va > vb
        return true
    end
    if va < vb
        return false
    end

    # if same column index => compare the density of rows
    va = length(a)
    vb = length(b)
    if va > vb
        return false
    end
    if va < vb
        return true
    end

    # hmm, equal rows?
    # should never happen
    return false
end

# Sort matrix upper rows (polynomial reducers) 
# by the leading column index and density.
#
# After the sort, the first row will have
# the left-most leading column index and, then, the smallest density.
function sort_matrix_upper_rows_decreasing!(matrix)
    #= smaller means pivot being more left  =#
    #= and density being smaller            =#

    inds = collect(1:(matrix.nup))
    cmp  = (x, y) -> matrix_row_decreasing_cmp(@inbounds(matrix.uprows[x]), @inbounds(matrix.uprows[y]))
    sort!(inds, lt=cmp, alg=_default_sorting_alg())

    matrix.uprows[1:(matrix.nup)] = matrix.uprows[inds]
    matrix.up2coef[1:(matrix.nup)] = matrix.up2coef[inds]

    matrix
end

# Sort matrix lower rows (polynomials to be reduced) 
# by the leading column index and density.
#
# After the sort, the first row will have
# the right-most leading column index and, then, the largest density.
function sort_matrix_lower_rows_increasing!(matrix)
    #= smaller means pivot being more right =#
    #= and density being larger             =#

    inds = collect(1:(matrix.nlow))
    cmp  = (x, y) -> matrix_row_increasing_cmp(@inbounds(matrix.lowrows[x]), @inbounds(matrix.lowrows[y]))
    sort!(inds, lt=cmp, alg=_default_sorting_alg())

    matrix.lowrows[1:(matrix.nlow)] = matrix.lowrows[inds]
    matrix.low2coef[1:(matrix.nlow)] = matrix.low2coef[inds]

    matrix
end

# Sorts the columns of the matrix (encoded in `col2hash` vector) 
# by the role of the corresponding monomial in the matrix,
# and then by the current monomial ordering decreasingly.
#
# See f4/matrix.jl for details
function sort_columns_by_hash!(col2hash::Vector{T}, symbol_ht::MonomialHashtable) where {T}
    hd = symbol_ht.hashdata
    es = symbol_ht.monoms

    function cmp(a, b, ord)
        @inbounds ha = hd[a]
        @inbounds hb = hd[b]
        if ha.idx != hb.idx
            return ha.idx > hb.idx
        end
        @inbounds ea = es[a]
        @inbounds eb = es[b]
        !monom_isless(ea, eb, ord)
    end

    ordcmp = (x, y) -> cmp(x, y, symbol_ht.ord)

    sort!(col2hash, lt=ordcmp, alg=_default_sorting_alg())
end

#------------------------------------------------------------------------------

# Given a vector of vectors of exponent vectors and coefficients,
# sort each vector wrt. the given monomial ordering `ord`.
# 
function sort_input_to_change_ordering!(
    exps::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    ord::AbstractMonomialOrdering
) where {M <: Monom, C <: Coeff}
    # for each polynomial encoded as 
    # two vectors exps[polyidx] and coeffs[polyidx]..
    @inbounds for polyidx in 1:length(exps)
        cmps =
            (x, y) ->
                monom_isless(@inbounds(exps[polyidx][y]), @inbounds(exps[polyidx][x]), ord)

        inds = collect(1:length(exps[polyidx]))

        sort!(inds, lt=cmps, alg=_default_sorting_alg())

        exps[polyidx][1:end] = exps[polyidx][inds]
        coeffs[polyidx][1:end] = coeffs[polyidx][inds]
    end
    nothing
end

#------------------------------------------------------------------------------

function sort_monom_indices_decreasing!(
    monoms::Vector{MonomIdx},
    cnt::Integer,
    ht::MonomialHashtable,
    ord::AbstractMonomialOrdering
)
    exps = ht.exponents

    cmps = (x, y) -> monom_isless(@inbounds(exps[y]), @inbounds(exps[x]), ord)

    sort_part!(monoms, 1, cnt, lt=cmps, alg=_default_sorting_alg())
end

function sort_term_indices_decreasing!(
    monoms::Vector{MonomIdx},
    coeffs::Vector{C},
    ht::MonomialHashtable,
    ord::AbstractMonomialOrdering
) where {C <: Coeff}
    exps = ht.exponents

    cmps =
        (x, y) -> monom_isless(@inbounds(exps[monoms[y]]), @inbounds(exps[monoms[x]]), ord)

    inds = collect(1:length(monoms))

    sort!(inds, lt=cmps, alg=Base.Sort.DEFAULT_UNSTABLE)

    monoms[1:end] = monoms[inds]
    coeffs[1:end] = coeffs[inds]
end
