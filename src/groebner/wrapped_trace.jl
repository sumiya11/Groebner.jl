# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# A wrapper around the F4 tracer exposed to the user

mutable struct WrappedTrace
    # One specialized trace per key.
    recorded_traces::Dict{Any, Any}   # Dict{Key, Trace}
    sys_support::Vector{Vector{Vector{IR_exponent}}}
    gb_support::Vector{Vector{Vector{IR_exponent}}}
    target_ord::AbstractMonomialOrdering

    WrappedTrace(sys_support, gb_support, target_ord) =
        new(Dict{Any, Any}(), sys_support, gb_support, target_ord)
end

function wrapped_trace_assert_compatible(trace::Trace, ring::PolyRing, key)
    @assert key <: CoeffGeneric ? trace.representation.coefftype <: CoeffGeneric :
            trace.representation.coefftype == key
    @assert trace.ring.ord == ring.ord
    @assert trace.ring.ground == ring.ground
    if key <: CoeffGeneric
        @assert typeof(trace.ring.characteristic) == typeof(ring.characteristic)
        @assert trace.ring.characteristic == ring.characteristic
    else
        @assert trace.ring.characteristic == typeof(trace.ring.characteristic)(ring.characteristic)
    end
    true
end

function wrapped_trace_save!(
    wrapped_trace::WrappedTrace,
    specialized_trace::Trace{CoeffType}
) where {CoeffType <: Coeff}
    key = CoeffType <: CoeffGeneric ? CoeffType : specialized_trace.representation.coefftype
    @assert key <: CoeffGeneric ? specialized_trace.representation.coefftype <: CoeffGeneric :
            specialized_trace.representation.coefftype == key
    if !isempty(wrapped_trace.recorded_traces)
        reference_trace = first(values(wrapped_trace.recorded_traces))
        @assert reference_trace.params.homogenize == specialized_trace.params.homogenize
        @assert reference_trace.ring.ord == specialized_trace.ring.ord
    end
    wrapped_trace.recorded_traces[key] = specialized_trace
end

function wrapped_trace_create_suitable_trace!(
    wrapped_trace::WrappedTrace,
    ring::PolyRing,
    params::AlgorithmParameters,
    ::Type{CoeffType}
) where {CoeffType <: Coeff}
    requested_key = CoeffType <: CoeffGeneric ? CoeffType : params.representation.coefftype
    # Try to find a suitable trace among the existing ones
    if haskey(wrapped_trace.recorded_traces, requested_key)
        trace = wrapped_trace.recorded_traces[requested_key]
        trace.ring.ord != ring.ord && throw(DomainError("ordering invalid in trace"))
        trace.ring = PolyRing(
            trace.ring.nvars,
            trace.ring.ord,
            typeof(trace.ring.characteristic)(ring.characteristic),
            trace.ring.ground
        )
        @invariant wrapped_trace_assert_compatible(trace, ring, requested_key)
        return trace
    end

    # Otherwise, create a new trace based on one of the existing ones
    default_trace = first(values(wrapped_trace.recorded_traces))
    default_trace.ring.ord != ring.ord && throw(DomainError("ordering invalid in trace"))
    if requested_key <: Union{CoeffGeneric, CoeffQQ}
        for trace in values(wrapped_trace.recorded_traces)
            if typeof(trace.ring.characteristic) == typeof(ring.characteristic)
                default_trace = trace
                break
            end
        end
    end
    new_trace = trace_copy(default_trace, params.representation, requested_key, ring.characteristic)
    new_char =
        requested_key <: Union{CoeffGeneric, CoeffQQ} ? ring.characteristic :
        convert(requested_key, ring.characteristic)
    new_trace.ring = PolyRing(ring.nvars, ring.ord, new_char, ring.ground)
    @invariant wrapped_trace_assert_compatible(new_trace, ring, requested_key)
    wrapped_trace_save!(wrapped_trace, new_trace)

    # the resulting trace may be in a invalid state, and needs to be filled with
    # the coefficients of the input polynomials
    new_trace
end

function wrapped_trace_check_input(
    trace::WrappedTrace,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}}
) where {I <: Integer, C <: Coeff}
    !(length(trace.sys_support) == length(monoms)) && return false
    for i in 1:length(monoms)
        !(length(trace.sys_support[i]) == length(monoms[i])) && return false
        for j in 1:length(monoms[i])
            if trace.sys_support[i][j] != monoms[i][j]
                return false
            end
        end
    end
    true
end

Base.show(io::IO, trace::WrappedTrace) = Base.show(io, MIME("text/plain"), trace)

function Base.show(io::IO, ::MIME"text/plain", wrapped_trace::WrappedTrace)
    println(io, "Recorded $(length(wrapped_trace.recorded_traces)) traces. Showing only one.\n")
    show(io, MIME("text/plain"), first(values(wrapped_trace.recorded_traces)))
end
