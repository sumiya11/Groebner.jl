# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backends for `groebner_learn` & `groebner_apply!`

###
# Learn stage

function groebner_learn0(polynomials, options)
    ring, monoms, coeffs =
        io_convert_polynomials_to_ir(polynomials, options)
    trace, gb_monoms, gb_coeffs = groebner_learn1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(
        ring, polynomials, gb_monoms, gb_coeffs, options)
    trace, result
end

# Proxy function for handling exceptions.
function groebner_learn1(ring::PolyRing, monoms, coeffs, options)
    try
        repr = io_select_polynomial_representation(ring, options)
        params = AlgorithmParameters(ring, repr, options)
        return _groebner_learn1(ring, monoms, coeffs, params, repr)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent."""
            repr = io_select_polynomial_representation(ring, options; large_exponents=true)
            params = AlgorithmParameters(ring, repr, options)
            return _groebner_learn1(ring, monoms, coeffs, params, repr)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function _groebner_learn1(ring, monoms, coeffs, params, repr)
    ring2, monoms2, coeffs2 = io_convert_ir_to_internal(ring, monoms, coeffs, params, repr)
    trace, gb_monoms2, gb_coeffs2 = groebner_learn2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)

    trace.representation = repr

    WrappedTraceF4(trace), gb_monoms, gb_coeffs
end

function groebner_learn2(ring, monoms, coeffs, params)
    if params.homogenize
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    trace, gb_monoms, gb_coeffs = _groebner_learn2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gbmonoms, gbcoeffs, params)
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

# Specialization for a single input
function groebner_apply0!(
    wrapped_trace::WrappedTraceF4,
    polynomials::AbstractVector,
    kws::KeywordArguments
)
    trace = get_trace!(wrapped_trace, polynomials, kws)
    @log :debug "Selected trace" trace.representation.coefftype

    flag, ring = io_extract_coeffs_raw!(trace, trace.representation, polynomials, kws)
    !flag && return (flag, polynomials)

    # TODO: this is a bit hacky
    params = AlgorithmParameters(
        ring,
        trace.representation,
        kws,
        orderings=(trace.params.original_ord, trace.params.target_ord)
    )
    ring = PolyRing(trace.ring.nvars, trace.ring.ord, ring.ch)

    flag, gb_monoms, gb_coeffs = groebner_apply1!(ring, trace, params)

    !flag && return (flag, polynomials)

    if trace.params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    flag, io_convert_to_output(ring, polynomials, gb_monoms, gb_coeffs, params)
end

# Specialization for a batch of several inputs
function groebner_apply0!(
    wrapped_trace::WrappedTraceF4,
    batch::NTuple{N, T},
    kws::KeywordArguments
) where {N, T <: AbstractVector}
    trace = get_trace!(wrapped_trace, batch, kws)
    @log :debug "Selected trace" trace.representation.coefftype

    flag, ring = io_io_extract_coeffs_raw_batched!(trace, trace.representation, batch, kws)
    !flag && return flag, batch

    # TODO: this is a bit hacky
    params = AlgorithmParameters(
        ring,
        trace.representation,
        kws,
        orderings=(trace.params.original_ord, trace.params.target_ord)
    )
    ring = PolyRing(trace.ring.nvars, trace.ring.ord, ring.ch)

    flag, gb_monoms, gb_coeffs = groebner_apply1!(ring, trace, params)

    !flag && return flag, batch

    if trace.params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    flag, io_convert_to_output_batched(ring, batch, gb_monoms, gb_coeffs, params)
end

function groebner_apply1!(ring, trace, params)
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

#=
Several assumptions are in place:
- input contains no zero polynomials (i.e., no empty vectors or vectors that
  contain only zeros),
- input coefficients are non-negative,
- input coefficients are smaller than the modulo.
=#
function groebner_applyX!(
    wrapped_trace::WrappedTraceF4,
    coeffs_zp::Vector{Vector{UInt32}},
    modulo::UInt32;
    options...
)
    kws = KeywordArguments(:groebner_apply!, options)

    logging_setup(kws)
    statistics_setup(kws)

    trace = get_trace!(wrapped_trace, modulo, kws)
    @log :debug "Selected trace" trace.representation.coefftype

    flag, ring = io_extract_coeffs_raw_X!(trace, trace.representation, coeffs_zp, modulo, kws)
    !flag && return (flag, coeffs_zp)

    # TODO: this is a bit hacky
    params = AlgorithmParameters(
        ring,
        trace.representation,
        kws,
        orderings=(trace.params.original_ord, trace.params.target_ord)
    )
    ring = PolyRing(trace.ring.nvars, trace.ring.ord, ring.ch)

    flag, gb_monoms, gb_coeffs = _groebner_apply1!(ring, trace, params)

    if trace.params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    flag, gb_coeffs
end
