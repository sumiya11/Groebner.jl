
mutable struct ComputationGraphF4{C, M, Ord}
    ring::PolyRing{Ord}
    input_basis::Basis{C}
    buf_basis::Basis{C}
    gb_basis::Basis{C}
    hashtable::MonomialHashtable{M, Ord}
    input_permutation::Vector{Int}
    # The number of columns,
    # The number of upper rows,
    # The number of lower rows
    matrix_infos::Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}
    matrix_nonzeroed_rows::Vector{Vector{Int}}
    matrix_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    matrix_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}
    # matrix_autored_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    # matrix_autored_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    # F4 iteration number --> index in the basis
    # matrix_tobereduced_rows::Vector{Vector{Int}}
    # matrix_tobereduced_mult::Vector{Vector{MonomIdx}}
    # # F4 iteration number --> index in the basis
    # matrix_reducers_rows::Vector{Vector{Int}}
    # matrix_reducers_mult::Vector{Vector{MonomIdx}}
    # matrix_columns::Vector{Vector{Int}}
end

function initialize_computation_graph_f4(ring, input_basis, gb_basis, hashtable, permutation)
    ComputationGraphF4(
        ring,
        input_basis,
        deepcopy_basis(gb_basis),
        gb_basis,
        hashtable,
        permutation,
        Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        # Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        # Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
    )
end

function finalize_graph!(graph::ComputationGraphF4) 
    # TODO: trim sizes
    graph.buf_basis = deepcopy_basis(graph.gb_basis)
    graph.buf_basis.ndivmasks = graph.input_basis.ndivmasks
    graph.buf_basis.nprocessed = graph.input_basis.nprocessed
    graph.buf_basis.ntotal = graph.input_basis.ntotal
    # TODO: set size to the allocated size
    # done.
    # graph.buf_basis.size = graph.input_basis.size
    nothing
end

function Base.show(io::IO, ::MIME"text/plain", graph::ComputationGraphF4)
    sz = Base.summarysize(graph) >> 10
    println(io, "Computation graph of F4 ($sz KiB)")
    println(io, "Recorded in TODO with parameters: TODO")
    total_iterations = length(graph.matrix_infos)
    total_matrix_low_rows = sum(x -> x.nlow, graph.matrix_infos)
    total_matrix_up_rows = sum(x -> x.nup, graph.matrix_infos)
    # TODO
    total_matrix_up_rows_useful = sum(length ∘ first, graph.matrix_upper_rows)
    total_matrix_low_rows_useful = sum(length ∘ first, graph.matrix_lower_rows)
    println(io, """
    Total iterations: $(total_iterations)
    Total upper matrix rows: $(total_matrix_up_rows) ($(round(total_matrix_up_rows_useful / total_matrix_up_rows * 100, digits=2)) % are useful)
    Total lower matrix rows: $(total_matrix_low_rows) ($(round(total_matrix_low_rows_useful / total_matrix_low_rows * 100, digits=2)) % are useful)""")
end
