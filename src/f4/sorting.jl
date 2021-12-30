

function sort_pairset_by_degree!(ps::Pairset, sz::Int)
    @info "Sorting by min degree"
    part = ps.pairs[1:sz]
    # TODO
    sort!(part, by=p -> p.deg)
    ps.pairs[1:sz] .= part
    ps
end

function sort_pairset_by_degree!(ps::Pairset)
    sort_pairset_by_degree!(ps, ps.load)
end

function sort_generators_by_something!(gens, load, basis)
    @info "Sorting first $load gens"
    part = gens[1:load]
    sort!(part, by=p -> length(basis.gens[p]))
    gens[1:load] .= part
    gens
end

function sort_columns_by_hash!(col2hash, symbol_ht)


    col2hash
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
    @info "Sorting input"
    part = basis.gens[1:basis.ntotal]
    sort!(part, lt=(x, y) -> exponent_isless(ht.exponents[x], ht.exponents[y]))
    basis.gens[1:basis.ntotal] .= part
    basis
end

#------------------------------------------------------------------------------

# by column index and density
function sort_matrix_rows_decreasing!(matrix)
    #= smaller means  pivot being more left  =#
    #= and density being smaller             =#


end

# by column index and density
function sort_matrix_rows_increasing!(matrix)
    #= smaller means  pivot being more right =#
    #= and density being larger              =#

end
