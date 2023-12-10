
# Information about F4 execution flow that can be later used for speeding up
# subsequent analogous computations.
mutable struct ComputationGraphF4{C, M, Ord, Ord2}
    time_spent::UInt64
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
    critical_pair_sequence::Vector{Tuple{Int, Int}}
    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}
    matrix_sorted_columns::Vector{Vector{Int}}
    params::AlgorithmParameters
    sweep_output::Bool
    representation::PolynomialRepresentation
    original_ord::Ord2
    term_sorting_permutations::Vector{Vector{Int}}
    term_homogenizing_permutations::Vector{Vector{Int}}
    homogenize::Bool

    input_signature::Vector{Int}
end

function initialize_computation_graph_f4(
    ring,
    input_basis,
    gb_basis,
    hashtable,
    permutation,
    params
)
    @log level = -5 "Initializing computation graph"
    input_signature = Vector{Int}(undef, input_basis.nfilled)
    @inbounds for i in 1:length(input_signature)
        input_signature[i] = length(input_basis.monoms[i])
    end
    ComputationGraphF4(
        time_ns(),
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
        Vector{Tuple{Int, Int}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Vector{Int}}(),
        params,
        params.sweep,
        PolynomialRepresentation(ExponentVector{UInt64}, UInt64),
        params.original_ord,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        params.homogenize,
        input_signature
    )
end

# TODO: is_applicable(graph, ideal...)

###
# Printing the graph

function finalize_graph!(graph::ComputationGraphF4)
    # TODO: trim sizes
    graph.buf_basis = deepcopy_basis(graph.gb_basis)
    graph.buf_basis.nnonredundant = graph.input_basis.nnonredundant
    graph.buf_basis.nprocessed = graph.input_basis.nprocessed
    graph.buf_basis.nfilled = graph.input_basis.nfilled
    graph.time_spent = time_ns() - graph.time_spent
    nothing
end

function Base.show(io::IO, ::MIME"text/plain", graph::ComputationGraphF4)
    sz = round((Base.summarysize(graph) / 2^20), digits=2)
    printstyled(
        io,
        "Computation graph ($sz MiB) recorded in $(round(graph.time_spent / 10^9, digits=3)) s.\n\n",
        bold=true
    )
    println(
        io,
        """
        Number of variables: $(graph.ring.nvars)
        Ground field characteristic: $(graph.ring.ch)
        Number of polynomials: $(graph.input_basis.nfilled) in input, $(graph.gb_basis.nfilled) in output
        """
    )

    permute_input =
        !isempty(graph.term_homogenizing_permutations) ||
        !isempty(graph.term_homogenizing_permutations)
    printstyled(io, "# Learn parameters\n", bold=true)
    println(
        io,
        """

        Original monomial ordering: $(graph.original_ord)
        Monom. representation: $(graph.representation.monomtype)
        Coeff. representation: $(graph.representation.coefftype)
        Sweep output: $(graph.sweep_output)
        Use homogenization: $(graph.homogenize)
        Permute input: $(permute_input)
        """
    )

    total_iterations = length(graph.matrix_infos)
    total_matrix_low_rows = sum(x -> x.nlow, graph.matrix_infos)
    total_matrix_up_rows = sum(x -> x.nup, graph.matrix_infos)
    total_matrix_up_rows_useful = sum(length ∘ first, graph.matrix_upper_rows)
    total_matrix_low_rows_useful = sum(length ∘ first, graph.matrix_lower_rows)
    critical_pair_degree_sequence = map(first, graph.critical_pair_sequence)
    critical_pair_count_sequence = map(last, graph.critical_pair_sequence)
    printstyled(io, "# F4 statistics\n", bold=true)
    println(
        io,
        """

        Iterations of F4: $(total_iterations)
        Monomial hashtable: $(graph.hashtable.load) / $(graph.hashtable.size) monomials filled
        Matrix total upper rows: $(total_matrix_up_rows) ($(round(total_matrix_up_rows_useful / total_matrix_up_rows * 100, digits=2)) % are useful)
        Matrix total lower rows: $(total_matrix_low_rows) ($(round(total_matrix_low_rows_useful / total_matrix_low_rows * 100, digits=2)) % are useful)
        Critical pair degree sequence: $(join(string.(critical_pair_degree_sequence), ","))
        Critical pair count sequence: $(join(string.(critical_pair_count_sequence), ","))"""
    )
end

function Base.show(io::IO, graph::ComputationGraphF4)
    Base.show(io, MIME("text/plain"), graph)
end
