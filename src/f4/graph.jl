
mutable struct ComputationGraphF4{C, M, Ord}
    ring::PolyRing{Ord}
    input_basis::Basis{C}
    buf_basis::Basis{C}
    gb_basis::Basis{C}
    hashtable::MonomialHashtable{M, Ord}
    input_permutation::Vector{Int}
    matrix_infos::Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}
    matrix_nonzeroed_rows::Vector{Vector{Int}}
    matrix_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    matrix_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}
    # hashtables_symbolic::Vector{MonomialHashtable{M, Ord}}
    matrix_sorted_columns::Vector{Vector{Int}}
end

function initialize_computation_graph_f4(
    ring,
    input_basis,
    gb_basis,
    hashtable,
    permutation
)
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
        Vector{Vector{Int}}()
        # Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        # Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
    )
end

function finalize_graph!(graph::ComputationGraphF4)
    # TODO: trim sizes
    graph.buf_basis = deepcopy_basis(graph.gb_basis)
    graph.buf_basis.nnonredundant = graph.input_basis.nnonredundant
    graph.buf_basis.nprocessed = graph.input_basis.nprocessed
    graph.buf_basis.nfilled = graph.input_basis.nfilled
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
    println(
        io,
        """
Total iterations: $(total_iterations)
Total upper matrix rows: $(total_matrix_up_rows) ($(round(total_matrix_up_rows_useful / total_matrix_up_rows * 100, digits=2)) % are useful)
Total lower matrix rows: $(total_matrix_low_rows) ($(round(total_matrix_low_rows_useful / total_matrix_low_rows * 100, digits=2)) % are useful)"""
    )
end
