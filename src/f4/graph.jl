
mutable struct ComputationGraphF4{C, M, Ord}
    input_basis::Basis{C}
    gb_basis::Basis{C}
    hashtable::MonomialHashtable{M, Ord}
    # The number of columns,
    # The number of upper rows,
    # The number of lower rows
    matrix_infos::Vector{Any}
    matrix_zeroed_rows::Vector{Vector{Int}}
    matrix_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    matrix_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    # F4 iteration number --> index in the basis
    # matrix_tobereduced_rows::Vector{Vector{Int}}
    # matrix_tobereduced_mult::Vector{Vector{MonomIdx}}
    # # F4 iteration number --> index in the basis
    # matrix_reducers_rows::Vector{Vector{Int}}
    # matrix_reducers_mult::Vector{Vector{MonomIdx}}
    # matrix_columns::Vector{Vector{Int}}
end

function initialize_computation_graph_f4(input_basis, gb_basis, hashtable)
    ComputationGraphF4(
        input_basis,
        gb_basis,
        hashtable,
        Vector{Any}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
    )
end

function Base.show(io::IO, ::MIME"text/plain", graph::ComputationGraphF4)
    println(io, "Computation graph of F4.")
    total_matrix_rows = sum(x -> x.nlow, graph.matrix_infos)
    total_redundant_rows = sum(length, graph.matrix_zeroed_rows)
    println(io, "Matrix rows in total: $(total_matrix_rows)")
    print("Redundant matrix rows: $(total_redundant_rows) ($(round(total_redundant_rows / total_matrix_rows * 100, digits=2)) %)")
end
