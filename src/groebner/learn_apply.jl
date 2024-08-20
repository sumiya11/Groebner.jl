# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backends for `groebner_learn` & `groebner_apply!`

###
# Learn stage

# polynomials => polynomials
function groebner_learn0(polynomials::AbstractVector, options::KeywordArguments)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    trace, gb_monoms, gb_coeffs = _groebner_learn1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    trace, result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function groebner_learn1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    _groebner_learn1(ring, monoms, coeffs, options)
end

function _groebner_learn1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    try
        params = AlgorithmParameters(ring, options)
        return __groebner_learn1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __groebner_learn1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __groebner_learn1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    term_sorting_permutations, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    trace, gb_monoms2, gb_coeffs2 = groebner_learn2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)

    trace.representation = params.representation
    trace.term_sorting_permutations = term_sorting_permutations

    WrappedTraceF4(trace), gb_monoms, gb_coeffs
end

# internal structs => internal structs
function groebner_learn2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return trace_initialize_empty(ring, monoms, coeffs, params),
        [monoms[1]],
        [coeffs[1]]
    end
    monoms, coeffs = _monoms, _coeffs

    term_homogenizing_permutation = Vector{Vector{Int}}()
    if params.homogenize
        term_homogenizing_permutation, ring, monoms, coeffs =
            homogenize_generators!(ring, monoms, coeffs, params)
    end

    trace, gb_monoms, gb_coeffs = _groebner_learn2(ring, monoms, coeffs, params)

    if params.homogenize
        trace.term_homogenizing_permutations = term_homogenizing_permutation
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    trace, gb_monoms, gb_coeffs
end

function _groebner_learn2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    trace, basis, pairset, hashtable =
        f4_initialize_structs_with_trace(ring, monoms, coeffs, params)
    f4_learn!(trace, ring, trace.gb_basis, pairset, hashtable, params)
    gb_monoms, gb_coeffs = basis_export_data(trace.gb_basis, trace.hashtable)
    trace, gb_monoms, gb_coeffs
end

###
# Apply stage

# polynomials => polynomials
function groebner_apply0!(
    wrapped_trace::WrappedTraceF4,
    polynomials::AbstractVector,
    options::KeywordArguments
)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    flag, gb_monoms, gb_coeffs =
        _groebner_apply1!(wrapped_trace, ring, monoms, coeffs, options)
    !flag && return (flag, polynomials)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    flag, result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function groebner_apply1!(
    wrapped_trace::WrappedTraceF4,
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    _groebner_apply1!(wrapped_trace, ring, monoms, coeffs, options)
end

function _groebner_apply1!(
    wrapped_trace::WrappedTraceF4,
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    trace = get_trace!(wrapped_trace, ring, options)

    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if trace.empty
        if isempty(_monoms)
            return true, [monoms[1]], [coeffs[1]]
        else
            return false, monoms, coeffs
        end
    end
    monoms, coeffs = _monoms, _coeffs

    # TODO: this is a bit hacky
    params = AlgorithmParameters(
        ring,
        options,
        orderings=(trace.params.original_ord, trace.params.target_ord)
    )
    params.representation = trace.representation
    ring = PolyRing(trace.ring.nvars, trace.ring.ord, ring.ch)

    flag = io_extract_coeffs_raw_X!(ring, trace, coeffs)
    !flag && return flag, [io_zero_monoms(Vector{UInt}, ring)], coeffs

    flag, gb_monoms2, gb_coeffs2 = groebner_apply2!(ring, trace, params)

    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring, gb_monoms2, gb_coeffs2, params)
    flag, gb_monoms, gb_coeffs
end

function groebner_apply2!(ring, trace, params)
    flag, gb_monoms, gb_coeffs = _groebner_apply2!(ring, trace, params)
    !flag && return (flag, gb_monoms, gb_coeffs)
    if trace.params.homogenize
        @assert false
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end
    flag, gb_monoms, gb_coeffs
end

function _groebner_apply2!(ring, trace, params)
    flag = f4_apply!(trace, ring, trace.buf_basis, params)

    gb_monoms, gb_coeffs = basis_export_data(trace.gb_basis, trace.hashtable)

    # Check once again that the sizes coincide
    if length(gb_monoms) != length(gb_coeffs)
        return false, gb_monoms, gb_coeffs
    end
    @inbounds for i in 1:length(gb_monoms)
        if length(gb_monoms[i]) != length(gb_coeffs[i])
            return false, gb_monoms, gb_coeffs
        end
    end

    flag, gb_monoms, gb_coeffs
end
