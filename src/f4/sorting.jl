# Sorting monomials, polynomials, and other things.

# Unstable algorithm should be a bit faster for large arrays
_default_sorting_alg() = Base.Sort.DEFAULT_UNSTABLE

# Sorts arr at the range of indices from..to 
function sort_part!(
    arr,
    from::Integer,
    to::Integer;
    lt=isless,
    alg=_default_sorting_alg(),
    by=identity
)
    # NOTE: this is perhaps type unstable
    ordr = Base.Sort.ord(lt, by, nothing)
    sort!(arr, from, to, alg, ordr)
    nothing
end

# Sorts polynomials from the basis by their leading monomial in the
# non-decreasing way by the given monomial ordering. Also sorts any arrays
# passed in the `abc` optional argument in the same order.
#
# Returns the sorting permutation.
function sort_polys_by_lead_increasing!(
    basis::Basis,
    hashtable::MonomialHashtable,
    abc...;
    ord::Ord=hashtable.ord
) where {Ord <: AbstractInternalOrdering}
    @log level = -4 "Sorting polynomials by their leading terms in non-decreasing order"
    b_monoms = basis.monoms
    h_monoms = hashtable.monoms
    permutation = collect(1:(basis.nfilled))
    cmps =
        (x, y) -> monom_isless(
            @inbounds(h_monoms[b_monoms[x][1]]),
            @inbounds(h_monoms[b_monoms[y][1]]),
            ord
        )
    # NOTE: stable sort to preserve the order of polynomials with the same lead
    sort!(permutation, lt=cmps, alg=Base.Sort.DEFAULT_STABLE)
    # use array assignment insted of elemewise assignment
    # (seems to compile to faster code)
    basis.monoms[1:(basis.nfilled)] = basis.monoms[permutation]
    basis.coeffs[1:(basis.nfilled)] = basis.coeffs[permutation]
    @inbounds for a in abc
        @invariant length(a) >= length(permutation)
        a[1:(basis.nfilled)] = a[permutation]
    end
    permutation
end

function is_sorted_by_lead_increasing(
    basis::Basis,
    hashtable::MonomialHashtable,
    ord::Ord=hashtable.ord
) where {Ord <: AbstractInternalOrdering}
    b_monoms = basis.monoms
    h_monoms = hashtable.monoms
    permutation = collect(1:(basis.nfilled))
    cmps =
        (x, y) -> monom_isless(
            @inbounds(h_monoms[b_monoms[x][1]]),
            @inbounds(h_monoms[b_monoms[y][1]]),
            ord
        )
    issorted(permutation, lt=cmps)
end

# Sorts critical pairs from the pairset in the range from..from+sz by the total
# degree of their lcms in increasing order
#
# Used in update once per one f4 iteration to sort pairs in pairset; Also used
# in normal selection strategy also once per iteration
function sort_pairset_by_degree!(ps::Pairset, from::Int, sz::Int)
    sort_part!(ps.pairs, from, from + sz, by=pair -> pair.deg)
end

# Sorts critical pairs from the pairset in the range from..from+sz by their
# sugar in increasing order
function sort_pairset_by_sugar!(
    ps::Pairset,
    from::Int,
    sz::Int,
    sugar_cubes::Vector{SugarCube}
)
    sort_part!(
        ps.pairs,
        from,
        from + sz,
        by=pair -> sugar_cubes[pair.poly1] + sugar_cubes[pair.poly2]
    )
    sugar = Vector{SugarCube}(undef, from + sz)
    @inbounds for i in 1:(from + sz)
        pair = ps.pairs[from + i - 1]
        sugar[i] = sugar_cubes[pair.poly1] + sugar_cubes[pair.poly2]
    end
    sugar
end

# Sorts the first `npairs` pairs from `pairset` in the non-decreasing order of
# their lcms by the given monomial ordering
function sort_pairset_by_lcm!(pairset::Pairset, npairs::Int, hashtable::MonomialHashtable)
    monoms = hashtable.monoms
    cmps =
        (pair1, pair2) -> monom_isless(
            @inbounds(monoms[pair1.lcm]),
            @inbounds(monoms[pair2.lcm]),
            hashtable.ord
        )
    sort_part!(pairset.pairs, 1, npairs, lt=cmps)
end

# Sorts polynomials by their position in the current basis (identity sort)
function sort_generators_by_position!(polys::Vector{Int}, load::Int)
    sort_part!(polys, 1, load)
end

###
# Sorting matrix rows and columns.
# See f4/matrix.jl for details.

# Compare sparse matrix rows a and b.
# A row is an array of integers, which are the indices of nonzero elements
function matrix_row_decreasing_cmp(a::Vector{T}, b::Vector{T}) where {T <: ColumnLabel}
    #= a, b - rows as arrays of nonzero indices =#
    # va and vb are the leading columns
    @inbounds va = a[1]
    @inbounds vb = b[1]
    if va > vb
        return false
    end
    if va < vb
        return true
    end
    # Unreachable.
    # If there are two rows in the upper part of the matrix with the same
    # leading term, something went wrong
    @invariant false
    # If the same leading column => compare the density of rows
    va = length(a)
    vb = length(b)
    if va > vb
        return true
    end
    if va < vb
        return false
    end
    false
end

# Compare sparse matrix rows a and b.
# A row is an array of integers, which are the indices of nonzero elements
function matrix_row_increasing_cmp(a::Vector{T}, b::Vector{T}) where {T <: ColumnLabel}
    #= a, b - rows as arrays of nonzero indices =#
    # va and vb are the leading columns
    @inbounds va = a[1]
    @inbounds vb = b[1]
    if va > vb
        return true
    end
    if va < vb
        return false
    end
    # If the same leading column => compare the density of rows
    va = length(a)
    vb = length(b)
    if va > vb
        return false
    end
    if va < vb
        return true
    end
    return false
end

# Sort matrix upper rows (polynomial reducers) by the leading column index and
# density.
#
# After the sort, the first (smallest) row will have the left-most leading
# column index and, then, the smallest density.
function sort_matrix_upper_rows!(matrix::MacaulayMatrix)
    #= smaller means pivot being more left  =#
    #= and density being smaller            =#
    permutation = collect(1:(matrix.nupper))
    cmp =
        (x, y) -> matrix_row_decreasing_cmp(
            @inbounds(matrix.upper_rows[x]),
            @inbounds(matrix.upper_rows[y])
        )
    sort!(permutation, lt=cmp, alg=_default_sorting_alg())
    matrix.upper_rows[1:(matrix.nupper)] = matrix.upper_rows[permutation]
    matrix.upper_to_coeffs[1:(matrix.nupper)] = matrix.upper_to_coeffs[permutation]
    # TODO: this is a bit hacky
    if !isempty(matrix.upper_to_mult)
        matrix.upper_to_mult[1:(matrix.nupper)] = matrix.upper_to_mult[permutation]
    end
    matrix
end

# Sort matrix lower rows (polynomials to be reduced) by the leading column index
# and density.
#
# After the sort, the first (smallest) row will have the right-most leading
# column index and, then, the largest density.
function sort_matrix_lower_rows!(matrix::MacaulayMatrix)
    #= smaller means pivot being more right =#
    #= and density being larger             =#
    permutation = collect(1:(matrix.nlower))
    cmp =
        (x, y) -> matrix_row_increasing_cmp(
            @inbounds(matrix.lower_rows[x]),
            @inbounds(matrix.lower_rows[y])
        )
    sort!(permutation, lt=cmp, alg=_default_sorting_alg())
    matrix.lower_rows[1:(matrix.nlower)] = matrix.lower_rows[permutation]
    matrix.lower_to_coeffs[1:(matrix.nlower)] = matrix.lower_to_coeffs[permutation]
    # TODO: this is a bit hacky
    if !isempty(matrix.lower_to_mult)
        matrix.lower_to_mult[1:(matrix.nlower)] = matrix.lower_to_mult[permutation]
    end
    matrix
end

# Sorts the columns of the matrix (encoded in `column_to_monom` vector) by the
# role of the corresponding monomial in the matrix, and then by the current
# monomial ordering decreasingly.
function sort_columns_by_labels!(
    column_to_monom::Vector{T},
    symbol_ht::MonomialHashtable
) where {T}
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

    sort!(column_to_monom, lt=ordcmp, alg=_default_sorting_alg())
end

# Given a vector of vectors of exponent vectors and coefficients, sort each
# vector wrt. the given monomial ordering `ord`.
#
# Returns the array of sorting permutations
function sort_input_terms_to_change_ordering!(
    exps::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    ord::AbstractInternalOrdering
) where {M <: Monom, C <: Coeff}
    permutations = Vector{Vector{Int}}(undef, length(exps))
    @inbounds for polyidx in 1:length(exps)
        comparator =
            (x, y) ->
                monom_isless(@inbounds(exps[polyidx][y]), @inbounds(exps[polyidx][x]), ord)

        permutation = collect(1:length(exps[polyidx]))
        sort!(permutation, lt=comparator, alg=_default_sorting_alg())

        exps[polyidx][1:end] = exps[polyidx][permutation]
        coeffs[polyidx][1:end] = coeffs[polyidx][permutation]

        permutations[polyidx] = permutation
    end
    permutations
end

function sort_monom_indices_decreasing!(
    monoms::Vector{MonomIdx},
    cnt::Integer,
    hashtable::MonomialHashtable,
    ord::AbstractInternalOrdering
)
    exps = hashtable.monoms

    cmps = (x, y) -> monom_isless(@inbounds(exps[y]), @inbounds(exps[x]), ord)

    sort_part!(monoms, 1, cnt, lt=cmps, alg=_default_sorting_alg())
end

function sort_term_indices_decreasing!(
    monoms::Vector{MonomIdx},
    coeffs::Vector{C},
    hashtable::MonomialHashtable,
    ord::AbstractInternalOrdering
) where {C <: Coeff}
    exps = hashtable.monoms

    cmps =
        (x, y) -> monom_isless(@inbounds(exps[monoms[y]]), @inbounds(exps[monoms[x]]), ord)

    inds = collect(1:length(monoms))

    sort!(inds, lt=cmps, alg=_default_sorting_alg())

    monoms[1:end] = monoms[inds]
    coeffs[1:end] = coeffs[inds]
end
