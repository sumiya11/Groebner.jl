


function sort_pairset_by_degree!(ps::Pairset)
    part = ps.pairs[1:ps.load]
    # TODO
    sort!(part, by=p -> p.deg)
    ps.pairs[1:ps.load] .= part
    ps
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
