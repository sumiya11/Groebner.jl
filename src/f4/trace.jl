# This file is a part of Groebner.jl. License is GNU GPL v2.

###
#  Tracing in f4

###
# The main struct

mutable struct TraceF4{C1 <: Coeff, C2 <: Coeff, M <: Monom, Ord1, Ord2}
    stopwatch_start::UInt64

    ring::PolyRing{Ord1, C2}
    original_ord::Ord2
    input_signature::Vector{Int}

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

    critical_pair_sequence::Vector{Tuple{Int, Int}}

    output_nonredundant_indices::Vector{Int}
    nonredundant_indices_before_reduce::Vector{Int}
    output_sort_indices::Vector{Int}

    params::AlgorithmParameters
    sweep_output::Bool
    representation::PolynomialRepresentation
    homogenize::Bool
end

function trace_initialize(
    ring::PolyRing,
    input_basis::Basis,
    gb_basis::Basis,
    hashtable::MonomialHashtable,
    permutation::Vector{Int},
    params::AlgorithmParameters
)
    @log :debug "Initializing the F4 trace"

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
        Vector{Tuple{Int, Int}}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        params,
        params.sweep,
        PolynomialRepresentation(ExponentVector{UInt64}, UInt64, false),
        params.homogenize
    )
end

function trace_deepcopy(
    trace::TraceF4{C1, C3, M, Ord1, Ord2}
) where {C1 <: Coeff, C3 <: Coeff, M <: Monom, Ord1, Ord2}
    # NOTE: does not provide the same guarantees as Base.deepcopy
    TraceF4(
        trace.stopwatch_start,
        PolyRing(trace.ring.nvars, trace.ring.ord, trace.ring.ch),
        deepcopy(trace.original_ord),
        copy(trace.input_signature),
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
        trace.homogenize
    )
end

function trace_copy(
    trace::TraceF4{C1, C3, M, Ord1, Ord2},
    ::Type{C2};
    deepcopy=false
) where {C1 <: Coeff, C3 <: Coeff, M <: Monom, Ord1, Ord2, C2 <: Coeff}
    new_sparse_row_coeffs = Vector{Vector{C2}}()
    new_input_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.input_basis, new_sparse_row_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.input_basis, new_sparse_row_coeffs)
    end

    new_buf_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.buf_basis.coeffs))
    for i in 1:length(trace.buf_basis.coeffs)
        !isassigned(trace.buf_basis.coeffs, i) && continue
        new_buf_basis_coeffs[i] = Vector{C2}(undef, length(trace.buf_basis.coeffs[i]))
    end
    new_buf_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.buf_basis, new_buf_basis_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.buf_basis, new_buf_basis_coeffs)
    end

    new_gb_basis_coeffs = Vector{Vector{C2}}(undef, length(trace.gb_basis.coeffs))
    for i in 1:length(trace.gb_basis.coeffs)
        !isassigned(trace.gb_basis.coeffs, i) && continue
        new_gb_basis_coeffs[i] = Vector{C2}(undef, length(trace.gb_basis.coeffs[i]))
    end
    new_gb_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(trace.gb_basis, new_gb_basis_coeffs)
    else
        basis_shallow_copy_with_new_coeffs(trace.gb_basis, new_gb_basis_coeffs)
    end

    new_representation = PolynomialRepresentation(trace.representation.monomtype, C2, false)
    new_ring = PolyRing(trace.ring.nvars, trace.ring.ord, zero(C2))

    TraceF4(
        trace.stopwatch_start,
        new_ring,
        trace.original_ord,
        trace.input_signature,
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
        trace.critical_pair_sequence,
        trace.output_nonredundant_indices,
        trace.nonredundant_indices_before_reduce,
        trace.output_sort_indices,
        trace.params,
        trace.sweep_output,
        new_representation,
        trace.homogenize
    )
end

function trace_finalize!(trace::TraceF4)
    # TODO: trim array sizes (Alex: well, this is whatever, memory usage won't be lowered)
    trace.buf_basis = basis_deepcopy(trace.gb_basis)
    trace.buf_basis.nnonredundant = trace.input_basis.nnonredundant
    trace.buf_basis.nprocessed = trace.input_basis.nprocessed
    trace.buf_basis.nfilled = trace.input_basis.nfilled
    trace.stopwatch_start = time_ns() - trace.stopwatch_start
    nothing
end

###
# A wrapper around the tracer exposed to the user

"""
    WrappedTraceF4


"""
mutable struct WrappedTraceF4
    recorded_traces::Dict{Any, Any}

    function WrappedTraceF4(
        trace::TraceF4{C1, C2, M, Ord1, Ord2}
    ) where {C1 <: Coeff, C2 <: Coeff, M <: Monom, Ord1, Ord2}
        recorded_traces = Dict{Any, Any}((C1, 42) => trace)
        WrappedTraceF4(recorded_traces)
    end

    function WrappedTraceF4(d::Dict{A, B}) where {A, B}
        new(d)
    end
end

"""
    trace_deepcopy(trace)

Returns a deep copy of `trace`. 
The original object and the copy are independent of each other.

Does not provide the same guarantees as `Base.deepcopy``.
"""
function trace_deepcopy(wrapped_trace::WrappedTraceF4)
    WrappedTraceF4(
        Dict(deepcopy(k) => trace_deepcopy(v) for (k, v) in wrapped_trace.recorded_traces)
    )
end

function get_default_trace(wrapped_trace::WrappedTraceF4)
    if length(wrapped_trace.recorded_traces) == 1
        # if there is only one trace stored
        first(values(wrapped_trace.recorded_traces))
    end

    for id in keys(wrapped_trace.recorded_traces)
        if id[1] <: Integer && id[2] == 42
            return wrapped_trace.recorded_traces[id]
        end
    end

    @unreachable
    first(values(wrapped_trace.recorded_traces))
end

function get_trace!(wrapped_trace::WrappedTraceF4, polynomials::AbstractVector, kwargs)
    ring = extract_ring(polynomials)
    get_trace!(wrapped_trace, ring.ch, kwargs)
end

function get_trace!(wrapped_trace::WrappedTraceF4, char, kwargs)
    trace = get_default_trace(wrapped_trace)

    # Fast path for the case when there exists a suitable trace
    coefftype = trace.representation.coefftype
    if (
        trace.representation.using_wide_type_for_coeffs && less_than_half(char, coefftype)
    ) || (!trace.representation.using_wide_type_for_coeffs && char <= typemax(coefftype))
        return trace
    end

    for id in keys(wrapped_trace.recorded_traces)
        if id[1] <: Integer && char <= typemax(id[1])
            return wrapped_trace.recorded_traces[id]
        end
    end

    # Handle the case when a wider coefficient type is required
    new_coefftype = io_get_tight_unsigned_int_type(char)
    @log :misc "Creating a new trace with coefficient type $new_coefftype"
    new_trace = trace_copy(trace, new_coefftype, deepcopy=false)
    wrapped_trace.recorded_traces[(new_coefftype, 0)] = new_trace

    new_trace
end

function get_trace!(
    wrapped_trace::WrappedTraceF4,
    batch::NTuple{N, T},
    kwargs
) where {N, T <: AbstractVector}
    # First, determine a suitable coefficient type for the given polynomials and
    # the learned trace
    default_trace = get_default_trace(wrapped_trace)
    coefftype = default_trace.representation.coefftype

    rings = map(extract_ring, batch)
    chars = map(ring -> ring.ch, rings)
    @log :misc """
    Determining a suitable coefficient type for the apply stage with characteristics $chars.
    On the learn stage, the coefficients were of type $coefftype, 
    and the $(typeof(default_trace.params.arithmetic)) arithmetic was used."""

    tight_signed_type = mapreduce(io_get_tight_signed_int_type, promote_type, chars)
    tight_unsigned_type = mapreduce(io_get_tight_unsigned_int_type, promote_type, chars)

    # The type of coefficients that will be used
    new_coefftype = if tight_signed_type == signed(tight_unsigned_type)
        tight_signed_type
    else
        @log :info """
        In the given batch of polynomials, the coefficient fields have
        characteristics $(chars) that do not fit into $(signed(tight_unsigned_type)),
        which may affect performance negatively.

        For best performance, use characteristics representable by $(Int32).
        Alternatively, please consider submitting a Github issue.
        """
        tight_unsigned_type
    end
    composite_coefftype = CompositeInt{N, new_coefftype}

    @log :misc """
    Will be storing polynomial coefficients as $composite_coefftype on the apply stage."""

    # Try to find a suitable trace among the existing ones
    for id in keys(wrapped_trace.recorded_traces)
        if id[1] <: composite_coefftype
            @log :misc "Re-using an existing trace with id = $id"
            return wrapped_trace.recorded_traces[id]
        end
    end

    # Otherwise, create a new trace based on one of the existing ones
    default_trace = get_default_trace(wrapped_trace)
    new_trace = trace_copy(default_trace, composite_coefftype, deepcopy=false)
    wrapped_trace.recorded_traces[(composite_coefftype, 1)] = new_trace

    # NOTE: the returned trace may be in a invalid state, and needs to be filled
    # with the coefficients of the input polynomials
    new_trace
end

###
# Printing the trace

function Base.show(io::IO, ::MIME"text/plain", wrapped_trace::WrappedTraceF4)
    println(
        io,
        "Recorded traces count: $(length(wrapped_trace.recorded_traces)). Printing the main one.\n"
    )
    show(io, MIME("text/plain"), get_default_trace(wrapped_trace))
end

function Base.show(io::IO, wrapped_trace::WrappedTraceF4)
    Base.show(io, MIME("text/plain"), wrapped_trace)
end

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
        Arithmetic type: $(typeof(trace.params.arithmetic))
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
