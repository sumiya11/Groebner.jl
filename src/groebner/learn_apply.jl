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
        ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    trace, gb_monoms2, gb_coeffs2 = groebner_learn2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = ir_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)

    trace.representation = params.representation
    trace.term_sorting_permutations = term_sorting_permutations
    trace.support = monoms

    WrappedTrace(trace), gb_monoms, gb_coeffs
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

    trace, gb_monoms, gb_coeffs = _groebner_learn2(ring, monoms, coeffs, params)

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
    wrapped_trace::WrappedTrace,
    polynomials::AbstractVector,
    options::KeywordArguments
)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    flag, gb_monoms, gb_coeffs =
        __groebner_apply1!(wrapped_trace, ring, monoms, coeffs, options)
    !flag && return (flag, polynomials)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    flag, result
end

# batch of polynomials => batch of polynomials
function groebner_apply_batch0!(
    wrapped_trace::WrappedTrace,
    batch::NTuple{N, T},
    options::KeywordArguments
) where {N, T}
    ir_batch = map(f -> io_convert_polynomials_to_ir(f, deepcopy(options)), batch)
    options = ir_batch[1][end]
    ir_batch = map(f -> f[1:3], ir_batch)
    flag, gb_batch = _groebner_apply_batch1!(wrapped_trace, ir_batch, options)
    !flag && return (flag, batch)
    result_ir = map(
        i -> io_convert_ir_to_polynomials(
            ir_batch[i][1],
            batch[i],
            gb_batch[i]...,
            options
        ),
        1:length(gb_batch)
    )
    flag, result_ir
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function groebner_apply1!(
    wrapped_trace::WrappedTrace,
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    __groebner_apply1!(wrapped_trace, ring, monoms, coeffs, options)
end

# batch of (exponent vectors, coefficients) 
# => 
# batch of (exponent vectors, coefficients)
function groebner_apply_batch1!(
    wrapped_trace::WrappedTrace,
    batch::NTuple{N, T},
    options::KeywordArguments
) where {N, T}
    batch = map(x -> ir_ensure_assumptions(x...), batch)
    _groebner_apply_batch1!(wrapped_trace, batch, options)
end

function _groebner_apply_batch1!(
    wrapped_trace::WrappedTrace,
    batch::NTuple{N, T},
    options::KeywordArguments
) where {N, T}
    flag, ring, monoms, coeffs = ir_pack_coeffs(batch)
    !flag && return flag, map(el -> el[2:end], batch)
    flag, gb_monoms, gb_coeffs =
        __groebner_apply1!(wrapped_trace, ring, monoms, coeffs, options)
    !flag && return flag, map(el -> el[2:end], batch)
    unpacked = ir_unpack_composite_coefficients(gb_coeffs)
    gb_batch = map(el -> (gb_monoms, el), unpacked)
    true, gb_batch
end

function __groebner_apply1!(
    wrapped_trace::WrappedTrace,
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    params = AlgorithmParameters(ring, options)

    trace = get_trace!(wrapped_trace, ring, params)

    flag = trace_check_input(trace, monoms, coeffs)
    !flag && return flag, monoms, coeffs

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

    flag = io_extract_coeffs_raw_X!(trace, coeffs)
    !flag && return flag, monoms, coeffs

    flag, gb_monoms2, gb_coeffs2 = groebner_apply2!(trace, params)

    gb_monoms, gb_coeffs = ir_convert_internal_to_ir(ring, gb_monoms2, gb_coeffs2, params)
    flag, gb_monoms, gb_coeffs
end

function groebner_apply2!(trace, params)
    flag, gb_monoms, gb_coeffs = _groebner_apply2!(trace, params)
    if !flag
        # Recover trace
        @log :info "Trace might be corrupted. Recovering..."
        trace.nfail += 1
        empty!(trace.matrix_sorted_columns)
        trace.buf_basis = basis_deepcopy(trace.input_basis)
        return flag, gb_monoms, gb_coeffs
    end
    flag, gb_monoms, gb_coeffs
end

function _groebner_apply2!(trace, params)
    flag = f4_apply!(trace, trace.ring, trace.buf_basis, params)

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
