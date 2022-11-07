
function monom_isless(a, b, ::Val{:degrevlex})
    exponent_isless_drl(a, b)
end

function monom_isless(a, b, ::Val{:lex})
    exponent_isless_lex(a, b)
end

function monom_isless(a, b, ::Val{:deglex})
    exponent_isless_dl(a, b)
end

#------------------------------------------------------------------------------

# sorts generators and corresponding coefficients from basis
# by their leading monomial in increasing ordering
#
# Used only once to sort input generators
function sort_gens_by_lead_increasing!(
            basis::Basis, ht::MonomialHashtable, abc...)
    gens = basis.monoms
    exps = ht.exponents

    inds = collect(1:basis.ntotal)

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[gens[x][1]]), 
        @inbounds(exps[gens[y][1]]), 
        Val(ht.ord)
    )

    sort!(inds, lt=cmp)

    # use array assignment insted of elemewise assignment
    basis.monoms[1:basis.ntotal]   = basis.monoms[inds]
    basis.coeffs[1:basis.ntotal] = basis.coeffs[inds]
    for a in abc
        a[1:basis.ntotal] = a[inds]
    end
end

function sort_gens_by_lead_increasing_in_reduce!(basis::Basis, ht::MonomialHashtable)
    gens = basis.monoms
    exps = ht.exponents
    nnr  = basis.nonred

    inds = collect(1:basis.nlead)

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[gens[nnr[x]][1]]), 
        @inbounds(exps[gens[nnr[y]][1]]), 
        Val(ht.ord)
    )
    sort!(inds, lt=cmp)

    # use array assignment insted of elemewise assignment
    basis.monoms[nnr[1:basis.nlead]]   = basis.monoms[nnr[inds]]
    basis.coeffs[nnr[1:basis.nlead]] = basis.coeffs[nnr[inds]]
end

function sort_gens_by_lead_increasing_in_standardize!(basis::Basis, ht::MonomialHashtable, ord)
    gens = basis.monoms
    exps = ht.exponents
    nnr  = basis.nonred

    inds = collect(1:basis.nlead)

    @assert inds == nnr

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[gens[x][1]]), 
        @inbounds(exps[gens[y][1]]), 
        Val(ord)
    )
    sort!(inds, lt=cmp)

    # use array assignment insted of elemewise assignment
    # Reason: slice array assignment uses fast setindex! based on unsafe copy
    basis.monoms[1:basis.nlead]   = basis.monoms[inds]
    basis.coeffs[1:basis.nlead] = basis.coeffs[inds]
    basis.lead[1:basis.nlead]   = basis.lead[inds]
end

#------------------------------------------------------------------------------

# Sorts pairs from pairset in range [from, from+size]
# by lcm total degrees in increasing order
#
# Used in update once per one f4 iteration to sort pairs in pairset
# Also used in normal selection strategy also once per iteration
function sort_pairset_by_degree!(ps::Pairset, from::Int, sz::Int)
    part = ps.pairs[from:from + sz]

    sort!(part, by=p -> p.deg)

    ps.pairs[from:from + sz] = part
end

#------------------------------------------------------------------------------

# sorts first npairs pairs from pairset by increasing of
# lcm exponent wrt the given monomial ordering
#
# Used only in normal selection strategy once per f4 iteration
function sort_pairset_by_lcm!(
        pairset::Pairset, npairs::Int, ht::MonomialHashtable)

    part = pairset.pairs[1:npairs]
    exps = ht.exponents

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[x.lcm]), 
        @inbounds(exps[y.lcm]), 
        Val(ht.ord)
    )
    sort!(part, lt=cmp)

    pairset.pairs[1:npairs] = part
end

#------------------------------------------------------------------------------

# sorts generators selected in normal strategy function
# by their ordering in the current basis (identity sort)
function sort_generators_by_position!(gens::Vector{Int}, load::Int)
    part = gens[1:load]

    sort!(part)

    gens[1:load] = part
end

#------------------------------------------------------------------------------

function matrix_row_decreasing_cmp(a, b)
    #= a, b - rows as arrays of exponent hashes =#

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
        return false
    end
    if va < vb
        return true
    end

    # hmm, equal rows?
    return false
end

function matrix_row_increasing_cmp(a, b)
    #= a, b - rows as arrays of exponent hashes =#

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
        return false
    end
    if va < vb
        return true
    end

    # hmm, equal rows?
    return false
end

# by column index and density
function sort_matrix_rows_decreasing!(matrix)
    #= smaller means  pivot being more left  =#
    #= and density being smaller             =#

    # @info "before up sort"
    #println(matrix.uprows)
    #println(matrix.up2coef)

    inds = collect(1:matrix.nup)
    cmp  = (x, y) ->  matrix_row_decreasing_cmp(
                                    @inbounds(matrix.uprows[x]),
                                    @inbounds(matrix.uprows[y]))
    sort!(inds, lt=cmp)

    matrix.uprows[1:matrix.nup]  = matrix.uprows[inds]
    matrix.up2coef[1:matrix.nup] = matrix.up2coef[inds]

    # @info "after up sort"
    #println(matrix.uprows)
    #println(matrix.up2coef)

    matrix
end

# by column index and density
function sort_matrix_rows_increasing!(matrix)
    #= smaller means  pivot being more right =#
    #= and density being larger =#

    # @info "before low sort"
    #println(matrix.lowrows)
    #println(matrix.low2coef)

    inds = collect(1:matrix.nlow)
    cmp  = (x, y) ->  matrix_row_increasing_cmp(
                                        @inbounds(matrix.lowrows[x]),
                                        @inbounds(matrix.lowrows[y]))
    sort!(inds, lt=cmp)

    matrix.lowrows[1:matrix.nlow]  = matrix.lowrows[inds]
    matrix.low2coef[1:matrix.nlow] = matrix.low2coef[inds]

    # @info "after low sort"
    #println(matrix.lowrows)
    #println(matrix.low2coef)

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

    sort!(col2hash, lt = ordcmp)
end

#------------------------------------------------------------------------------
# TODO: ordering??

function sort_input_to_change_ordering(exps, coeffs, ord::Symbol)
    for polyidx in 1:length(exps)
        cmp = (x, y) -> monom_isless(
            @inbounds(exps[polyidx][y]), 
            @inbounds(exps[polyidx][x]), 
            Val(ord)
        )

        inds = collect(1:length(exps[polyidx]))

        sort!(inds, lt=cmp)

        exps[polyidx][1:end] = exps[polyidx][inds]
        coeffs[polyidx][1:end] = coeffs[polyidx][inds]
    end
end

#------------------------------------------------------------------------------

function sort_monoms_increasing!(monoms::Vector{MonomIdx}, cnt, ht, ord::Symbol)
    exps = ht.exponents

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[monoms[x]]), 
        @inbounds(exps[monoms[y]]), 
        Val(ht.ord)
    )

    inds = collect(1:cnt)

    sort!(inds, lt=cmp)

    monoms[1:cnt] = monoms[inds]
end

function sort_monoms_decreasing!(monoms::Vector{MonomIdx}, cnt, ht, ord::Symbol)
    exps = ht.exponents

    cmp = (x, y) -> monom_isless(
        @inbounds(exps[monoms[y]]), 
        @inbounds(exps[monoms[x]]), 
        Val(ord)
    )

    inds = collect(1:cnt)

    sort!(inds, lt=cmp)

    monoms[1:cnt] = monoms[inds]
end

function sort_terms_decreasing!(monoms::Vector{MonomIdx}, coeffs, ht, ord::Symbol)
    exps = ht.exponents
    
    cmp = (x, y) -> monom_isless(
        @inbounds(exps[monoms[y]]), 
        @inbounds(exps[monoms[x]]), 
        Val(ord)
    )

    inds = collect(1:length(monoms))

    sort!(inds, lt=cmp)

    monoms[1:end] = monoms[inds]
    coeffs[1:end] = coeffs[inds]
end
