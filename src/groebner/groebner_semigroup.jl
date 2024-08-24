
function groebner_semigroup(polynomials, varmap, relations, groups, options)
    isempty(polynomials) && throw(DomainError("Empty input."))
    # @log :info "" polynomials relations groups
    _options = deepcopy(options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    ring2, monoms2, coeffs2, options2 =
        io_convert_polynomials_to_ir(relations, deepcopy(_options))
    # TODO: multihom
    ((ring3, monoms3, coeffs3, _), (ring4, monoms4, coeffs4, _)) =
        map(group -> io_convert_polynomials_to_ir(group, deepcopy(_options)), groups)
    ring5, monoms5, coeffs5, options5 =
        io_convert_polynomials_to_ir(varmap, deepcopy(_options))
    params = AlgorithmParameters(ring, options)
    params.linalg = LinearAlgebra(:deterministic, :sparse)
    params.reduced = false
    params.representation = PolynomialRepresentation(
        ExponentVector{UInt32},
        params.representation.coefftype,
        params.representation.using_wide_type_for_coeffs
    )
    # @log :info "" ring ring2 ring3 ring4 ring5
    @assert ring == ring2 == ring3 == ring4 == ring5
    gb_monoms, gb_coeffs = __groebner_semigroup1(
        ring,
        monoms,
        coeffs,
        params,
        monoms2,
        coeffs2,
        options2,
        monoms3,
        coeffs3,
        monoms4,
        coeffs4,
        monoms5,
        coeffs5
    )
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    result
end

function __groebner_semigroup1(
    ring,
    monoms,
    coeffs,
    params::AlgorithmParameters,
    monoms2,
    coeffs2,
    options2,
    monoms3,
    coeffs3,
    monoms4,
    coeffs4,
    monoms5,
    coeffs5
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring, monoms, coeffs =
        ir_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    _, _, monoms2, coeffs2 =
        ir_convert_ir_to_internal(ring, monoms2, coeffs2, params, params.representation)
    _, _, monoms3, coeffs3 =
        ir_convert_ir_to_internal(ring, monoms3, coeffs3, params, params.representation)
    _, _, monoms4, coeffs4 =
        ir_convert_ir_to_internal(ring, monoms4, coeffs4, params, params.representation)
    _, _, monoms5, coeffs5 =
        ir_convert_ir_to_internal(ring, monoms5, coeffs5, params, params.representation)
    SEMIGROUP_RELATIONS[] = monoms2
    SEMIGROUP_GROUPS[] = [monoms3, monoms4]
    SEMIGROUP_VARMAP[] = monoms5
    # @log 100 "" monoms monoms2
    # @log 100 "" SEMIGROUP_RELATIONS[]
    # @log 100 "" SEMIGROUP_GROUPS[]
    # @log 100 "" SEMIGROUP_VARMAP[]
    gb_monoms, gb_coeffs =
        _groebner_semigroup2(ring, vcat(monoms, monoms2), vcat(coeffs, coeffs2), params)
    gb_monoms, gb_coeffs = ir_convert_internal_to_ir(ring, gb_monoms, gb_coeffs, params)
    gb_monoms, gb_coeffs
end

function _groebner_semigroup2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    f4!(ring, basis, pairset, hashtable, params)
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    gbmonoms, gbcoeffs
end
