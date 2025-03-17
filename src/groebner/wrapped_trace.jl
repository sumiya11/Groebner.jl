# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# A wrapper around the F4 tracer exposed to the user

mutable struct WrappedTrace
    # One specialized trace per key.
    recorded_traces::Dict{Any, Any}   # Dict{Key, Trace}
    sys_support::Vector{Vector{Vector{IR_exponent}}}
    gb_support::Vector{Vector{Vector{IR_exponent}}}

    WrappedTrace(sys_support, gb_support) = new(Dict{Any, Any}(), sys_support, gb_support)
end

function wrapped_trace_save!(wrapped_trace::WrappedTrace, specialized_trace)
    wrapped_trace.recorded_traces[specialized_trace.representation.coefftype] = specialized_trace
end

function wrapped_trace_create_suitable_trace!(
    wrapped_trace::WrappedTrace,
    ring::PolyRing,
    params::AlgorithmParameters
)
    # Try to find a suitable trace among the existing ones
    if haskey(wrapped_trace.recorded_traces, params.representation.coefftype)
        trace = wrapped_trace.recorded_traces[params.representation.coefftype]
        trace.ring = PolyRing(
            trace.ring.nvars,
            trace.ring.ord,
            typeof(trace.ring.characteristic)(ring.characteristic)
        )
        return trace
    end

    # Otherwise, create a new trace based on one of the existing ones
    default_trace = first(values(wrapped_trace.recorded_traces))
    new_trace = trace_copy(default_trace, params.representation)
    new_trace.ring.ord != ring.ord && throw(DomainError("ordering invalid in trace"))
    new_trace.ring = PolyRing(
        ring.nvars,
        ring.ord,
        convert(params.representation.coefftype, ring.characteristic),
        ring.ground
    )
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
