# This file is a part of Groebner.jl. License is GNU GPL v2.

# Tracer for F4

mutable struct Trace{C1 <: Coeff, C2 <: Coeff, M <: Monom, Ord1, Ord2}
    stopwatch::UInt64
    empty::Bool

    ring::PolyRing{Ord1, C2}
    original_ord::Ord2
    support::Vector{Vector{Vector{Int}}}

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
    matrix_sorted_columns::Vector{Vector{Int}}
    matrix_pivot_signatures::Vector{UInt64}
    matrix_pivot_indices::Vector{Vector{Int}}
    matrix_is_columns_cached::Bool

    critical_pair_sequence::Vector{Tuple{Int, Int}}

    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}

    params::AlgorithmParameters
    sweep_output::Bool
    representation::PolynomialRepresentation
    homogenize::Bool

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
    hashtable = hashtable_initialize(ring, params.rng, M, 2)
    trace = trace_initialize(ring, basis, basis, hashtable, Int[], params)
    trace.empty = true
    trace
end

function trace_initialize(
    ring::PolyRing,
    input_basis::Basis,
    gb_basis::Basis,
    hashtable::MonomialHashtable,
    permutation::Vector{Int},
    params::AlgorithmParameters
)
    Trace(
        time_ns(),
        false,
        ring,
        params.original_ord,
        Vector{Vector{Vector{Int}}}(),
        input_basis,
        basis_deepcopy(gb_basis),
        gb_basis,
        hashtable,
        permutation,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        Vector{NamedTuple{(:nup, :nlow, :ncols), Tuple{Int64, Int64, Int64}}}(),
        Vector{Vector{Int}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomId}}}(),
        Vector{Tuple{Vector{Int}, Vector{MonomId}}}(),
        Vector{Vector{Int}}(),
        Vector{UInt64}(),
        Vector{Vector{Int}}(),
        false,
        Vector{Tuple{Int, Int}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        params,
        params.sweep,
        PolynomialRepresentation(ExponentVector{UInt64}, UInt64, false),
        params.homogenize,
        0,
        0
    )
end

function trace_deepcopy(
    trace::Trace{C1, C3, M, Ord1, Ord2}
) where {C1 <: Coeff, C3 <: Coeff, M <: Monom, Ord1, Ord2}
    # NOTE: does not provide the same guarantees as Base.deepcopy
    Trace(
        trace.stopwatch,
        trace.empty,
        PolyRing(trace.ring.nvars, trace.ring.ord, trace.ring.ch),
        deepcopy(trace.original_ord),
        deepcopy(trace.support),
        basis_deepcopy(trace.input_basis),
        basis_deepcopy(trace.buf_basis),
        basis_deepcopy(trace.gb_basis),
        # we are assuming that the hashtable is frozen at this point!
        # maybe set a flag in the hashtable?
        deepcopy(trace.hashtable),
        copy(trace.input_permutation),
        map(copy, trace.term_sorting_permutations),
        map(copy, trace.term_homogenizing_permutations),
        # matrix information must not be mutated no matter what
        trace.matrix_infos,
        trace.matrix_nonzeroed_rows,
        trace.matrix_upper_rows,
        trace.matrix_lower_rows,
        map(copy, trace.matrix_sorted_columns),
        copy(trace.matrix_pivot_signatures),
        map(copy, trace.matrix_pivot_indices),
        trace.matrix_is_columns_cached,
        copy(trace.critical_pair_sequence),
        copy(trace.output_nonredundant_indices),
        copy(trace.nonredundant_indices_before_reduce),
        copy(trace.output_sort_indices),
        deepcopy(trace.params),
        trace.sweep_output,
        PolynomialRepresentation(
            trace.representation.monomtype,
            trace.representation.coefftype,
            trace.representation.using_wide_type_for_coeffs
        ),
        trace.homogenize,
        trace.napply,
        trace.nfail
    )
end

function trace_copy(
    trace::Trace{C1, C3, M, Ord1, Ord2},
    ::Type{C2},
    using_wide_type_for_coeffs::Bool;
    deepcopy=false
) where {C1 <: Coeff, C3 <: Coeff, M <: Monom, Ord1, Ord2, C2 <: Coeff}
    new_sparse_row_coeffs = Vector{Vector{C2}}()
    new_input_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.input_basis, new_sparse_row_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.input_basis, new_sparse_row_coeffs)
    end

    new_buf_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.buf_basis.coeffs))
    # for i in 1:length(trace.buf_basis.coeffs)
    for i in 1:(trace.buf_basis.nfilled)
        # !isassigned(trace.buf_basis.coeffs, i) && continue
        new_buf_basis_coeffs[i] = Vector{C2}(undef, length(trace.buf_basis.coeffs[i]))
    end
    new_buf_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.buf_basis, new_buf_basis_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.buf_basis, new_buf_basis_coeffs)
    end

    new_gb_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.gb_basis.coeffs))
    # for i in 1:length(trace.gb_basis.coeffs)
    for i in 1:(trace.gb_basis.nfilled)
        # !isassigned(trace.gb_basis.coeffs, i) && continue
        new_gb_basis_coeffs[i] = Vector{C2}(undef, length(trace.gb_basis.coeffs[i]))
    end
    new_gb_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.gb_basis, new_gb_basis_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.gb_basis, new_gb_basis_coeffs)
    end

    new_representation = PolynomialRepresentation(
        trace.representation.monomtype,
        C2,
        using_wide_type_for_coeffs
    )
    new_ring = PolyRing(trace.ring.nvars, trace.ring.ord, zero(C2))

    Trace(
        trace.stopwatch,
        trace.empty,
        new_ring,
        trace.original_ord,
        trace.support,
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
        trace.sweep_output,
        new_representation,
        trace.homogenize,
        trace.napply,
        trace.nfail
    )
end

function trace_finalize!(trace::Trace)
    trace.buf_basis = basis_deepcopy(trace.gb_basis)
    trace.buf_basis.nnonredundant = trace.input_basis.nnonredundant
    trace.buf_basis.nprocessed = trace.input_basis.nprocessed
    trace.buf_basis.nfilled = trace.input_basis.nfilled
    trace.stopwatch = time_ns() - trace.stopwatch
    nothing
end

function trace_check_input(
    trace::Trace,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}}
) where {I <: Integer, C <: Coeff}
    !(length(trace.support) == length(monoms)) && return false
    for i in 1:length(monoms)
        !(length(trace.support[i]) == length(monoms[i])) && return false
        for j in 1:length(monoms[i])
            if trace.support[i][j] != monoms[i][j]
                return false
            end
        end
    end
    true
end

###
# A wrapper around the tracer exposed to the user

mutable struct WrappedTrace
    # For each type of coefficients, we maintain a separate tracer object
    recorded_traces::Dict{Any, Any}

    function WrappedTrace(
        trace::Trace{C1, C2, M, Ord1, Ord2}
    ) where {C1 <: Coeff, C2 <: Coeff, M <: Monom, Ord1, Ord2}
        recorded_traces = Dict{Any, Any}((C1, 42) => trace)
        WrappedTrace(recorded_traces)
    end

    function WrappedTrace(d::Dict{A, B}) where {A, B}
        new(d)
    end
end

# Does not provide the same guarantees as Base.deepcopy.
function trace_deepcopy(wrapped_trace::WrappedTrace)
    WrappedTrace(
        Dict(deepcopy(k) => trace_deepcopy(v) for (k, v) in wrapped_trace.recorded_traces)
    )
end

function get_default_trace(wrapped_trace::WrappedTrace)
    if length(wrapped_trace.recorded_traces) == 1
        # if there is only one trace stored
        first(values(wrapped_trace.recorded_traces))
    end

    for id in keys(wrapped_trace.recorded_traces)
        if id[1] <: Integer && id[2] == 42
            return wrapped_trace.recorded_traces[id]
        end
    end

    first(values(wrapped_trace.recorded_traces))
end

function get_trace!(
    wrapped_trace::WrappedTrace,
    ring::PolyRing,
    params::AlgorithmParameters
)
    # Try to find a suitable trace among the existing ones
    for id in keys(wrapped_trace.recorded_traces)
        if id[1] == params.representation.coefftype
            @log :misc "Re-using an existing trace with id = $id"
            trace = wrapped_trace.recorded_traces[id]
            trace.ring.ch = ring.ch
            return trace
        end
    end

    # Otherwise, create a new trace based on one of the existing ones
    default_trace = get_default_trace(wrapped_trace)
    new_trace = trace_copy(
        default_trace,
        params.representation.coefftype,
        params.representation.using_wide_type_for_coeffs,
        deepcopy=false
    )
    new_trace.ring =
        PolyRing(ring.nvars, ring.ord, convert(params.representation.coefftype, ring.ch))
    wrapped_trace.recorded_traces[(params.representation.coefftype, 1)] = new_trace

    # the resulting trace may be in a invalid state, and needs to be filled with
    # the coefficients of the input polynomials
    new_trace
end

###
# Printing the trace

Base.show(io::IO, trace::WrappedTrace) = Base.show(io, MIME("text/plain"), trace)
Base.show(io::IO, trace::Trace) = Base.show(io, MIME("text/plain"), trace)

function Base.show(io::IO, ::MIME"text/plain", wrapped_trace::WrappedTrace)
    println(
        io,
        """Recorded $(length(wrapped_trace.recorded_traces)) traces with IDs: $(collect(keys(wrapped_trace.recorded_traces)))
        Showing only one.\n"""
    )
    show(io, MIME("text/plain"), get_default_trace(wrapped_trace))
end

function Base.show(io::IO, ::MIME"text/plain", trace::Trace)
    tm = round(trace.stopwatch / 10^9, digits=3)
    sz = round((Base.summarysize(trace) / 2^20), digits=2)
    printstyled(io, "# Trace of F4 recorded in $(tm) s ($sz MiB).\n", bold=true)
    println(
        io,
        """
        ring  : Z[x1,...,x$(trace.ring.nvars)] mod $(trace.ring.ch)
        input : $(trace.input_basis.nfilled) polynomials
        output: $(trace.gb_basis.nfilled) polynomials
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
        sweep        : $(trace.sweep_output)
        homogenize   : $(trace.homogenize)
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
        matrix largest : $((0, 0))
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
