# This file is a part of Groebner.jl. License is GNU GPL v2.

# Tracer for F4

mutable struct Trace{C1 <: Coeff, C2 <: Coeff, M <: Monom, Ord1, Ord2}
    stopwatch::UInt64
    empty::Bool

    ring::PolyRing{Ord1, C2}
    original_ord::Ord2

    # Buffers for storing basis elements
    input_basis::Basis{C1}
    buf_basis::Basis{C1}
    gb_basis::Basis{C1}
    hashtable::MonomialHashtable{M, Ord1}

    # Permutation of input polynomials and their terms, if needed
    input_permutation::Vector{Int}
    term_sorting_permutations::Vector{Vector{Int}}
    term_homogenizing_permutations::Vector{Vector{Int}}

    # Information about the non-redundant rows of the matrix
    matrix_infos::Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}
    matrix_nonzeroed_rows::Vector{Vector{Int}}
    matrix_upper_rows::Vector{Tuple{Vector{Int}, Vector{MonomId}}}
    matrix_lower_rows::Vector{Tuple{Vector{Int}, Vector{MonomId}}}
    matrix_sorted_columns::Vector{Vector{Int32}}
    matrix_pivot_signatures::Vector{UInt64}
    matrix_pivot_indices::Vector{Vector{Int}}
    matrix_is_columns_cached::Bool

    critical_pair_sequence::Vector{Tuple{Int, Int}}

    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}

    params::AlgorithmParameters
    representation::PolynomialRepresentation

    napply::Int
    nfail::Int
end

function trace_initialize_empty(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis = basis_initialize(ring, 0, C)
    hashtable = hashtable_initialize(ring, params.rng, M; use_divmask=params.use_divmask)
    trace = trace_initialize(ring, basis, hashtable, Int[], params)
    trace.empty = true
    trace
end

function trace_initialize(
    ring::PolyRing,
    input_basis::Basis,
    hashtable::MonomialHashtable,
    permutation::Vector{Int},
    params::AlgorithmParameters
)
    Trace(
        time_ns(),
        false,
        ring,
        params.original_ord,
        basis_deepcopy(input_basis),
        basis_deepcopy(input_basis),
        input_basis,
        hashtable,
        permutation,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomId}}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomId}}}(),
        Vector{Vector{Int32}}(),
        Vector{UInt64}(),
        Vector{Vector{Int}}(),
        false,
        Vector{Tuple{Int, Int}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        params,
        PolynomialRepresentation(ExponentVector{UInt64}, UInt64, false),
        0,
        0
    )
end

function trace_copy(
    trace::Trace{C1, C3, M, Ord1, Ord2},
    repr::PolynomialRepresentation
) where {C1 <: Coeff, C3 <: Coeff, M <: Monom, Ord1, Ord2}
    C2 = repr.coefftype
    new_sparse_row_coeffs = Vector{Vector{C2}}()
    new_input_basis = basis_shallow_copy_with_new_coeffs(trace.input_basis, new_sparse_row_coeffs)

    new_buf_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.buf_basis.coeffs))
    # for i in 1:length(trace.buf_basis.coeffs)
    for i in 1:(trace.buf_basis.n_filled)
        # !isassigned(trace.buf_basis.coeffs, i) && continue
        new_buf_basis_coeffs[i] = Vector{C2}(undef, length(trace.buf_basis.coeffs[i]))
    end
    new_buf_basis = basis_shallow_copy_with_new_coeffs(trace.buf_basis, new_buf_basis_coeffs)

    new_gb_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.gb_basis.coeffs))
    # for i in 1:length(trace.gb_basis.coeffs)
    for i in 1:(trace.gb_basis.n_filled)
        # !isassigned(trace.gb_basis.coeffs, i) && continue
        new_gb_basis_coeffs[i] = Vector{C2}(undef, length(trace.gb_basis.coeffs[i]))
    end
    new_gb_basis = basis_shallow_copy_with_new_coeffs(trace.gb_basis, new_gb_basis_coeffs)

    new_representation = PolynomialRepresentation(
        trace.representation.monomtype,
        C2,
        repr.using_wide_type_for_coeffs
    )
    new_ring = PolyRing(trace.ring.nvars, trace.ring.ord, zero(C2), trace.ring.ground)

    Trace(
        trace.stopwatch,
        trace.empty,
        new_ring,
        trace.original_ord,
        new_input_basis,
        new_buf_basis,
        new_gb_basis,
        trace.hashtable,
        trace.input_permutation,
        trace.term_sorting_permutations,
        trace.term_homogenizing_permutations,
        trace.matrix_infos,
        trace.matrix_nonzeroed_rows,
        trace.matrix_upper_rows,
        trace.matrix_lower_rows,
        copy(trace.matrix_sorted_columns),
        trace.matrix_pivot_signatures,
        trace.matrix_pivot_indices,
        trace.matrix_is_columns_cached,
        trace.critical_pair_sequence,
        trace.output_nonredundant_indices,
        trace.nonredundant_indices_before_reduce,
        trace.output_sort_indices,
        trace.params,
        new_representation,
        trace.napply,
        trace.nfail
    )
end

function trace_finalize!(trace::Trace)
    trace.buf_basis = basis_deepcopy(trace.gb_basis)
    trace.buf_basis.n_nonredundant = trace.input_basis.n_nonredundant
    trace.buf_basis.n_processed = trace.input_basis.n_processed
    trace.buf_basis.n_filled = trace.input_basis.n_filled
    trace.stopwatch = time_ns() - trace.stopwatch
    nothing
end

function trace_export_matrices(trace::Trace)
    @assert length(trace.matrix_infos) ==
            length(trace.matrix_lower_rows) ==
            length(trace.matrix_upper_rows)
    matrices = []
    nvars = trace.ring.nvars
    for i in 1:length(trace.matrix_infos)
        reducers = (
            index=trace.matrix_upper_rows[i][1],
            multiplier=map(
                m -> monom_to_vector!(zeros(Int, nvars), trace.hashtable.monoms[m]),
                trace.matrix_upper_rows[i][2]
            )
        )
        to_be_reduced = (
            index=trace.matrix_lower_rows[i][1],
            multiplier=map(
                m -> monom_to_vector!(zeros(Int, nvars), trace.hashtable.monoms[m]),
                trace.matrix_lower_rows[i][2]
            )
        )
        matrix = (reducers=reducers, to_be_reduced=to_be_reduced)
        push!(matrices, matrix)
    end
    # if the last iteration is interreduction
    if length(trace.critical_pair_sequence) + 1 == length(matrices)
        matrices[end] = (reducers=matrices[end].to_be_reduced, to_be_reduced=matrices[end].reducers)
    end
    map(identity, matrices)
end

Base.show(io::IO, trace::Trace) = Base.show(io, MIME("text/plain"), trace)

function Base.show(io::IO, ::MIME"text/plain", trace::Trace)
    tm = round(trace.stopwatch / 10^9, digits=3)
    sz = round((Base.summarysize(trace) / 2^20), digits=2)
    printstyled(io, "# Trace of F4 recorded in $(tm) s ($sz MiB).\n", bold=true)
    println(
        io,
        """
        ring  : Z[x1,...,x$(trace.ring.nvars)] mod $(trace.ring.characteristic)
        input : $(trace.input_basis.n_filled) polynomials
        output: $(trace.gb_basis.n_filled) polynomials
        apply : $(trace.napply) / $(trace.nfail) (success/fail)
        """
    )

    permute_input =
        !isempty(trace.term_homogenizing_permutations) ||
        !isempty(trace.term_homogenizing_permutations)
    printstyled(io, "# Parameters\n", bold=true)
    println(
        io,
        """
        input order  : $(trace.original_ord)
        output order : $(trace.ring.ord)
        homogenize   : $(trace.params.homogenize)
        permute      : $(permute_input)
        monom. type  : $(trace.representation.monomtype)
        coeff. type  : $(trace.representation.coefftype)
        arithmetic   : $(typeof(trace.params.arithmetic))
        """
    )

    total_iterations = length(trace.matrix_infos)
    total_matrix_low_rows = sum(x -> x.nlow, trace.matrix_infos; init=0)
    total_matrix_up_rows = sum(x -> x.nup, trace.matrix_infos; init=0)
    total_matrix_up_rows_useful = sum(length ∘ first, trace.matrix_upper_rows; init=0)
    total_matrix_low_rows_useful = sum(length ∘ first, trace.matrix_lower_rows; init=0)
    critical_pair_degree_sequence = map(first, trace.critical_pair_sequence)
    critical_pair_count_sequence = map(last, trace.critical_pair_sequence)
    printstyled(io, "# F4 statistics\n", bold=true)
    print(
        io,
        """
        iterations     : $(total_iterations)
        hashtable      : $(trace.hashtable.load) / $(trace.hashtable.size) filled
        matrix up-rows : $(total_matrix_up_rows) ($(round(total_matrix_up_rows_useful / total_matrix_up_rows * 100, digits=2)) % useful)
        matrix low-rows: $(total_matrix_low_rows) ($(round(total_matrix_low_rows_useful / total_matrix_low_rows * 100, digits=2)) % useful)
        pair degrees   : """
    )
    if length(critical_pair_degree_sequence) > 0
        print(io, string(critical_pair_degree_sequence[1]))
        print(io, ",")
    end
    for i in 2:length(critical_pair_degree_sequence)
        if critical_pair_degree_sequence[i] < critical_pair_degree_sequence[i - 1]
            printstyled(io, string(critical_pair_degree_sequence[i]), color=:red, bold=true)
        else
            print(io, string(critical_pair_degree_sequence[i]))
        end
        if i != length(critical_pair_degree_sequence)
            print(io, ",")
        end
    end
    println(io)
    println(
        io,
        """
        pair count     : $(join(string.(critical_pair_count_sequence), ","))"""
    )
end
