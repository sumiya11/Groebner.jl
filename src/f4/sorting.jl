
_default_alg() = Base.Sort.DEFAULT_UNSTABLE

# Sorts vector v at indices from:to
function sort_part!(v, from::Integer, to::Integer;
                    lt=isless, 
                    alg=_default_alg(), 
                    by=identity,
                    rev=false)
    ordr = Base.Sort.ord(lt,by,rev,Base.Sort.Forward)
    sort!(v, from, to, alg, ordr)
end

function choose_comparator(ord::Symbol)
    if ord === :degrevlex
        exponent_isless_drl
    elseif ord === :deglex
        exponent_isless_dl
    elseif ord === :lex
        exponent_isless_lex
    end
end

#------------------------------------------------------------------------------

# sorts polynomials from the basis
# by their leading monomial in the non-decreasing way by term ordering
function sort_gens_by_lead_increasing!(
            basis::Basis, ht::MonomialHashtable, abc...; ord::Symbol=ht.ord)
    gens = basis.monoms
    exps = ht.exponents

    inds = collect(1:basis.ntotal)

    cmp = choose_comparator(ord)
    cmps = (x, y) -> cmp(
        @inbounds(exps[gens[x][1]]), 
        @inbounds(exps[gens[y][1]])
    )

    sort!(inds, lt=cmps, alg=_default_alg())

    # use array assignment insted of elemewise assignment
    basis.monoms[1:basis.ntotal] = basis.monoms[inds]
    basis.coeffs[1:basis.ntotal] = basis.coeffs[inds]
    for a in abc
        a[1:basis.ntotal] = a[inds]
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
    sort_part!(
        ps.pairs, 
        from, from+sz,
        by=p -> p.deg, 
        alg=_default_alg()
    )
end

#------------------------------------------------------------------------------

# sorts first `npairs` pairs from `pairset` by non-decreasing of
# lcm exponent vector wrt the given monomial ordering
function sort_pairset_by_lcm!(
        pairset::Pairset, npairs::Int, ht::MonomialHashtable)

    exps = ht.exponents

    cmp = choose_comparator(ht.ord)
    cmps = (x, y) -> cmp(
        @inbounds(exps[x.lcm]), 
        @inbounds(exps[y.lcm])
    )

    sort_part!(
        pairset.pairs, 
        1, npairs,
        lt=cmps, 
        alg=_default_alg()
    )
end

#------------------------------------------------------------------------------

# sorts generators selected in normal strategy function
# by their ordering in the current basis (identity sort)
function sort_generators_by_position!(gens::Vector{Int}, load::Int)
    sort_part!(
        gens, 
        1, load,
        alg=_default_alg()
    )
end

#------------------------------------------------------------------------------

function matrix_row_decreasing_cmp(a, b)
    #= a, b - rows as arrays of exponent indices =#
    @inbounds va = a[1]
    @inbounds vb = b[1]

    if va > vb
        return false
    end
    if va < vb
        return true
    end

    # same column index => compare density of rows
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

function matrix_row_increasing_cmp(a, b)
    #= a, b - rows as arrays of exponent hashes =#
    @inbounds va = a[1]
    @inbounds vb = b[1]

    if va > vb
        return true
    end
    if va < vb
        return false
    end

    # same column index => compare density of rows
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

# by column index and density
function sort_matrix_rows_decreasing!(matrix)
    #= smaller means  pivot being more left  =#
    #= and density being smaller             =#

    inds = collect(1:matrix.nup)
    cmp  = (x, y) ->  matrix_row_decreasing_cmp(
                                    @inbounds(matrix.uprows[x]),
                                    @inbounds(matrix.uprows[y]))
    sort!(inds, lt=cmp, alg=_default_alg())

    matrix.uprows[1:matrix.nup]  = matrix.uprows[inds]
    matrix.up2coef[1:matrix.nup] = matrix.up2coef[inds]

    matrix
end

# by column index and density
function sort_matrix_rows_increasing!(matrix)
    #= smaller means  pivot being more right =#
    #= and density being larger =#

    inds = collect(1:matrix.nlow)
    cmp  = (x, y) ->  matrix_row_increasing_cmp(
                                        @inbounds(matrix.lowrows[x]),
                                        @inbounds(matrix.lowrows[y]))
    sort!(inds, lt=cmp, alg=_default_alg())

    matrix.lowrows[1:matrix.nlow]  = matrix.lowrows[inds]
    matrix.low2coef[1:matrix.nlow] = matrix.low2coef[inds]

    matrix
end

function sort_columns_by_hash!(col2hash, symbol_ht)
    hd = symbol_ht.hashdata
    es = symbol_ht.exponents

    function cmp(a, b, exponent_isless)
        @inbounds ha = hd[a]
        @inbounds hb = hd[b]
        if ha.idx != hb.idx
            return ha.idx > hb.idx
        end
        @inbounds ea = es[a]
        @inbounds eb = es[b]
        !exponent_isless(ea, eb)
    end

    if symbol_ht.ord == :degrevlex
        ordcmp = (x, y) -> cmp(x, y, exponent_isless_drl)
    elseif symbol_ht.ord == :lex
        ordcmp = (x, y) -> cmp(x, y, exponent_isless_lex)
    else
        ordcmp = (x, y) -> cmp(x, y, exponent_isless_dl)
    end

    sort!(col2hash, lt = ordcmp, alg=_default_alg())
end

#------------------------------------------------------------------------------

function sort_input_to_change_ordering(exps, coeffs, ord::Symbol)
    for polyidx in 1:length(exps)
        cmp = choose_comparator(ord)
        cmps = (x, y) -> cmp(
            @inbounds(exps[polyidx][y]), 
            @inbounds(exps[polyidx][x]), 
        )

        inds = collect(1:length(exps[polyidx]))

        sort!(inds, lt=cmps, alg=_default_alg())

        exps[polyidx][1:end] = exps[polyidx][inds]
        coeffs[polyidx][1:end] = coeffs[polyidx][inds]
    end
end

#------------------------------------------------------------------------------

function sort_monoms_decreasing!(monoms::Vector{MonomIdx}, cnt, ht, ord::Symbol)
    exps = ht.exponents

    cmp = choose_comparator(ord)
    cmps = (x, y) -> cmp(
        @inbounds(exps[y]), 
        @inbounds(exps[x])
    )

    sort_part!(
        monoms, 
        1, cnt,
        lt=cmps,
        alg=_default_alg()
    )
end

function sort_terms_decreasing!(monoms::Vector{MonomIdx}, coeffs, ht, ord::Symbol)
    exps = ht.exponents
    
    cmp = choose_comparator(ord)
    cmps = (x, y) -> cmp(
        @inbounds(exps[monoms[y]]), 
        @inbounds(exps[monoms[x]])
    )

    inds = collect(1:length(monoms))

    sort!(inds, lt=cmps, alg=Base.Sort.DEFAULT_UNSTABLE)

    monoms[1:end] = monoms[inds]
    coeffs[1:end] = coeffs[inds]
end
