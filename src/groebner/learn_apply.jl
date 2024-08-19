# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backends for `groebner_learn` & `groebner_apply!`

###
# Learn stage

function groebner_learn0(polynomials, options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    trace, gb_monoms, gb_coeffs = groebner_learn1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    trace, result
end

# Proxy function for handling exceptions.
function groebner_learn1(ring::PolyRing, monoms, coeffs, options)
    try
        params = AlgorithmParameters(ring, options)
        return _groebner_learn1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent."""
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return _groebner_learn1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function _groebner_learn1(ring, monoms, coeffs, params)
    @invariant ir_is_valid(ring, monoms, coeffs)
    term_sorting_permutations, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    trace, gb_monoms2, gb_coeffs2 = groebner_learn2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)

    trace.representation = params.representation
    trace.term_sorting_permutations = term_sorting_permutations

    WrappedTraceF4(trace), gb_monoms, gb_coeffs
end

function groebner_learn2(ring, monoms, coeffs, params)
    monoms = filter(!isempty, monoms)
    coeffs = filter(!isempty, coeffs)

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
    ring,
    monoms,
    coeffs::Vector{Vector{C}},
    params
) where {C <: CoeffZp}
    @log :misc "Groebner learn phase over Z_p"
    # Initialize F4 structs
    trace, basis, pairset, hashtable =
        f4_initialize_structs_with_trace(ring, monoms, coeffs, params)
    @log :all "Before F4:" basis
    f4_learn!(trace, ring, trace.gb_basis, pairset, hashtable, params)
    @log :all "After F4:" basis
    gb_monoms, gb_coeffs = basis_export_data(trace.gb_basis, trace.hashtable)
    trace, gb_monoms, gb_coeffs
end

###
# Apply stage

function groebner_apply0!(
    wrapped_trace::WrappedTraceF4,
    polynomials::AbstractVector,
    options::KeywordArguments
)
    # flag, ring, coeffs = io_extract_coeffs_raw!(trace, trace.representation, polynomials, options)
    # ring, coeffs = io_extract_coeffs_ir!(trace, trace.representation, polynomials, options)
    # TODO: too slow
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    flag, gb_monoms, gb_coeffs = groebner_apply1!(wrapped_trace, ring, coeffs, options)
    !flag && return (flag, polynomials)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    flag, result
end

function groebner_apply1!(wrapped_trace, ring, coeffs, options)
    trace = get_trace!(wrapped_trace, ring, options)
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
    @log :misc "Groebner Apply phase"
    @log :misc "Applying modulo $(ring.ch)"

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
