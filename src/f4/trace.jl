# Tracing for f4

mutable struct TraceF4{C <: Coeff, M <: Monom, Ord1, Ord2}
    stopwatch_start::UInt64

    ring::PolyRing{Ord1}
    original_ord::Ord2
    input_signature::Vector{Int}

    input_basis::Basis{C}
    buf_basis::Basis{C}
    gb_basis::Basis{C}

    hashtable::MonomialHashtable{M, Ord1}

    # Permutation of input polynomials and their terms, if needed
    input_permutation::Vector{Int}
    term_sorting_permutations::Vector{Vector{Int}}
    term_homogenizing_permutations::Vector{Vector{Int}}

    # Information about the non-redundant rows of the matrix
    matrix_infos::Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}
    matrix_nonzeroed_rows::Vector{Vector{Int}}
    matrix_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    matrix_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}
    matrix_sorted_columns::Vector{Vector{Int}}

    critical_pair_sequence::Vector{Tuple{Int, Int}}

    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}

    params::AlgorithmParameters
    sweep_output::Bool
    representation::PolynomialRepresentation
    homogenize::Bool
end

function initialize_trace_f4(
    ring::PolyRing,
    input_basis::Basis,
    gb_basis::Basis,
    hashtable::MonomialHashtable,
    permutation::Vector{Int},
    params::AlgorithmParameters
)
    @log level = -5 "Initializing the F4 trace"

    input_signature = Vector{Int}(undef, input_basis.nfilled)
    @inbounds for i in 1:length(input_signature)
        input_signature[i] = length(input_basis.monoms[i])
    end

    TraceF4(
        time_ns(),
        ring,
        params.original_ord,
        input_signature,
        input_basis,
        deepcopy_basis(gb_basis),
        gb_basis,
        hashtable,
        permutation,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomIdx}}}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Int, Int}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        params,
        params.sweep,
        PolynomialRepresentation(ExponentVector{UInt64}, UInt64),
        params.homogenize
    )
end

function finalize_trace!(trace::TraceF4)
    # TODO: trim array sizes
    trace.buf_basis = deepcopy_basis(trace.gb_basis)
    trace.buf_basis.nnonredundant = trace.input_basis.nnonredundant
    trace.buf_basis.nprocessed = trace.input_basis.nprocessed
    trace.buf_basis.nfilled = trace.input_basis.nfilled
    trace.stopwatch_start = time_ns() - trace.stopwatch_start
    nothing
end

###
# Printing the trace

function Base.show(io::IO, ::MIME"text/plain", trace::TraceF4)
    sz = round((Base.summarysize(trace) / 2^20), digits=2)
    printstyled(
        io,
        "Trace of F4 ($sz MiB) recorded in $(round(trace.stopwatch_start / 10^9, digits=3)) s.\n\n",
        bold=true
    )
    println(
        io,
        """
        Number of variables: $(trace.ring.nvars)
        Ground field characteristic: $(trace.ring.ch)
        Number of polynomials: $(trace.input_basis.nfilled) in input, $(trace.gb_basis.nfilled) in output
        """
    )

    permute_input =
        !isempty(trace.term_homogenizing_permutations) ||
        !isempty(trace.term_homogenizing_permutations)
    printstyled(io, "# Learn parameters\n", bold=true)
    println(
        io,
        """

        Original monomial ordering: $(trace.original_ord)
        Monom. representation: $(trace.representation.monomtype)
        Coeff. representation: $(trace.representation.coefftype)
        Sweep output: $(trace.sweep_output)
        Use homogenization: $(trace.homogenize)
        Permute input: $(permute_input)
        """
    )

    total_iterations = length(trace.matrix_infos)
    total_matrix_low_rows = sum(x -> x.nlow, trace.matrix_infos)
    total_matrix_up_rows = sum(x -> x.nup, trace.matrix_infos)
    total_matrix_up_rows_useful = sum(length ∘ first, trace.matrix_upper_rows)
    total_matrix_low_rows_useful = sum(length ∘ first, trace.matrix_lower_rows)
    critical_pair_degree_sequence = map(first, trace.critical_pair_sequence)
    critical_pair_count_sequence = map(last, trace.critical_pair_sequence)
    printstyled(io, "# F4 statistics\n", bold=true)
    println(
        io,
        """

        Iterations of F4: $(total_iterations)
        Monomial hashtable: $(trace.hashtable.load) / $(trace.hashtable.size) monomials filled
        Matrix total upper rows: $(total_matrix_up_rows) ($(round(total_matrix_up_rows_useful / total_matrix_up_rows * 100, digits=2)) % are useful)
        Matrix total lower rows: $(total_matrix_low_rows) ($(round(total_matrix_low_rows_useful / total_matrix_low_rows * 100, digits=2)) % are useful)
        Critical pair degree sequence: $(join(string.(critical_pair_degree_sequence), ","))
        Critical pair count sequence: $(join(string.(critical_pair_count_sequence), ","))"""
    )
end

function Base.show(io::IO, trace::TraceF4)
    Base.show(io, MIME("text/plain"), trace)
end
