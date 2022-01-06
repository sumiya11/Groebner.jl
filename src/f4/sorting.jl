

function sort_pairset_by_degree!(ps::Pairset, from::Int, sz::Int)
    @info "Sorting by min degree"
    part = ps.pairs[from:sz]
    # TODO
    sort!(part, by=p -> p.deg)
    ps.pairs[from:sz] .= part
    ps
end

function sort_pairset_by_degree!(ps::Pairset)
    sort_pairset_by_degree!(ps, 1, ps.load)
end

function sort_generators_by_position!(gens, load)
    @info "Sorting first $load gens"
    part = gens[1:load]
    sort!(part)
    gens[1:load] .= part
    gens
end

function sort_pairset_by_lcm!(pairset, npairs, ht)
    part = pairset.pairs[1:npairs]
    sort!(part, lt=(x, y) -> exponent_isless(ht.exponents[x.lcm], ht.exponents[y.lcm]))
    pairset.pairs[1:npairs] .= part
    pairset
end

function sort_columns_by_hash!(col2hash, symbol_ht)
    function cmp(x, y, ht)
        ha = ht.hashdata[x]
        hb = ht.hashdata[y]

        exponent_isless(ht.exponents[x], ht.exponents[y])
    end

    closedcmp = (x, y) -> cmp(x, y, symbol_ht)
    sort!(col2hash, lt=closedcmp)
end

#------------------------------------------------------------------------------

#=
    :degrevlex exponent vector comparison
=#
function exponent_isless(e1, e2)
    n = length(e1)
    for i in 1:n
        if e1[i] > e2[i]
            return false
        end
    end
    return e1[end] < e2[end]
end

function sort_gens_by_lead_increasing!(basis, ht)
    @info "Sorting input" basis.gens basis.coeffs

    gens = basis.gens

    inds = collect(1:basis.ntotal)
    cmp  = (x, y) -> exponent_isless(ht.exponents[gens[x][1]], ht.exponents[gens[y][1]])
    sort!(inds, lt=cmp)

    basis.gens[1:basis.ntotal] .= basis.gens[1:basis.ntotal][inds]
    basis.coeffs[1:basis.ntotal] .= basis.coeffs[1:basis.ntotal][inds]

    @info "Sorted input" basis.gens basis.coeffs

    basis
end

#------------------------------------------------------------------------------

function matrix_row_decreasing_cmp(a, b)
    #= a, b - rows as arrays of exponent hashes =#

    va = a[1]
    vb = b[1]

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
        return true
    end
    if va < vb
        return false
    end

    # hmm, equal rows?
    return false
end

function matrix_row_increasing_cmp(a, b)
    #= a, b - rows as arrays of exponent hashes =#

    va = a[1]
    vb = b[1]

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

    @info "before up sort"
    println(matrix.uprows)
    println(matrix.up2coef)

    inds = collect(1:matrix.nup)
    cmp  = (x, y) ->  matrix_row_decreasing_cmp(matrix.uprows[x], matrix.uprows[y])
    sort!(inds, lt=cmp)

    matrix.uprows[1:matrix.nup] .= matrix.uprows[1:matrix.nup][inds]
    matrix.up2coef[1:matrix.nup] .= matrix.up2coef[1:matrix.nup][inds]

    @info "after up sort"
    println(matrix.uprows)
    println(matrix.up2coef)

    matrix
end

# by column index and density
function sort_matrix_rows_increasing!(matrix)
    #= smaller means  pivot being more right =#
    #= and density being larger =#

    @info "before low sort"
    println(matrix.lowrows)
    println(matrix.low2coef)

    inds = collect(1:matrix.nlow)
    cmp  = (x, y) ->  matrix_row_increasing_cmp(matrix.lowrows[x], matrix.lowrows[y])
    sort!(inds, lt=cmp)

    matrix.lowrows[1:matrix.nlow] .= matrix.lowrows[1:matrix.nlow][inds]
    matrix.low2coef[1:matrix.nlow] .= matrix.low2coef[1:matrix.nlow][inds]

    @info "after low sort"
    println(matrix.lowrows)
    println(matrix.low2coef)

    matrix
end

function sort_gens_terms_decreasing!(basis, ht, i)
    inds = collect(1:length(basis.gens[i]))
    gen = basis.gens[i]
    cmp  = (x, y) -> exponent_isless(ht.exponents[gen[x]], ht.exponents[gen[y]])
    sort!(inds, lt=cmp)
    basis.gens[i] = basis.gens[i][inds]
    basis.coeffs[i] = basis.coeffs[i][inds]
    basis
end
