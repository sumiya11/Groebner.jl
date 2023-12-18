# Input-output conversions of polynomials.

# Our conventions on some edge cases:
# - Trying to compute a Groebner basis of an empty set is an error.
# - The Groebner basis of [0] is [0].
# - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]

# NOTE: Polynomials from frontend packages, such as AbstractAlgebra.jl and
# DynamicPolynomials.jl, may not support some of the orderings supported by
# Groebner.jl. We compute the basis in the requested ordering in Groebner.jl,
# and then order the terms in the output according to the ordering of the input
# polynomials

# NOTE: internally, 0 polynomial is represented with an empty vector of
# monomials and an empty vector of coefficients

@noinline function __throw_inexact_coeff_conversion(c, T)
    throw(
        DomainError(
            c,
            """
            Coefficient $c in the output basis cannot be converted exactly to $T. 
            Using big arithmetic in the input could fix this."""
        )
    )
end

@noinline function __throw_input_not_supported(msg, val)
    throw(DomainError(val, "This type of input is not supported, sorry.\n$msg"))
end

"""
    PolyRing

Polynomial ring of computation.
"""
mutable struct PolyRing{
    Ord <: Union{AbstractMonomialOrdering, AbstractInternalOrdering},
    C <: Union{CoeffZp, CompositeCoeffZp}
}
    # Number of variables
    nvars::Int
    # Monomial ordering
    ord::Ord
    # Characteristic of the coefficient ring
    ch::C
end

###
# Selecting polynomial representation

"""
    PolynomialRepresentation

Internal representation of polynomials.
"""
struct PolynomialRepresentation
    monomtype::Type
    coefftype::Type
    # NOTE: If this field is false, then any implementation of the arithmetic in
    # Z/Zp must cast the coefficients into a wider integer type before
    # performing any arithmetic operations to avoid the risk of overflow.
    using_wide_type_for_coeffs::Bool
end

"""
    select_polynomial_representation(polynomials, keywords; hint=:none)

Given an array of input polynomials tries to select a suitable representation
for coefficients and exponent vectors.

Additionally, `hint` can be specified to one of the following:

- `:large_exponents`: use at least 32 bits per exponent.
"""
function select_polynomial_representation(
    polynomials::AbstractVector,
    kws::KeywordsHandler;
    hint::Symbol=:none
)
    if !(hint in (:none, :large_exponents))
        @log level = 1_000 "The given hint=$hint was discarded"
    end
    frontend, npolys, char, nvars, ordering = peek_at_polynomials(polynomials)
    monomtype = select_monomtype(char, npolys, nvars, ordering, kws, hint)
    coefftype, using_wide_type_for_coeffs =
        select_coefftype(char, npolys, nvars, ordering, kws, hint)
    basering = iszero(char) ? :qq : :zp
    @log level = -1 "Frontend: $frontend"
    @log level = -1 """
    Input: $npolys polynomials over $basering in $nvars variables
    Ordering: $ordering"""
    @log level = -1 """
    Internal representation: 
    monomials are $monomtype
    coefficients are $coefftype
    wide type for coefficients: $using_wide_type_for_coeffs"""
    PolynomialRepresentation(monomtype, coefftype, using_wide_type_for_coeffs)
end

function select_monomtype(char, npolys, nvars, ordering, kws, hint)
    @log level = -1 """
    Selecting monomial representation.
    Given hint hint=$hint. Keyword argument monoms=$(kws.monoms)"""
    if hint === :large_exponents
        @log level = -1 "As hint=$hint was provided, using 64 bits per single exponent"
        # use 64 bits if large exponents detected
        desired_monom_type = ExponentVector{UInt64}
        @assert is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
    end

    # If homogenization is requested, or if a part of the ordering is
    # lexicographical, the ideal will potentially be homogenized later.
    if kws.homogenize === :yes || (
        kws.homogenize === :auto && (
            kws.ordering isa Lex ||
            kws.ordering isa ProductOrdering ||
            (ordering == :lex && kws.ordering isa InputOrdering)
        )
    )
        @log level = -1 """
        As homogenize=:yes/:auto was provided, 
        representing monomials as dense exponent vectors"""
        desired_monom_type = ExponentVector{UInt32}
        @assert is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
    end

    ExponentSize = UInt8
    variables_per_word = div(sizeof(UInt), sizeof(ExponentSize))
    # if dense representation is requested
    if kws.monoms === :dense
        @assert is_supported_ordering(ExponentVector{ExponentSize}, kws.ordering)
        return ExponentVector{ExponentSize}
    end
    # if sparse representation is requested
    if kws.monoms === :sparse
        if is_supported_ordering(
            SparseExponentVector{ExponentSize, nvars, Int},
            kws.ordering
        )
            return SparseExponentVector{ExponentSize, nvars, Int}
        end
        @log level = 1 """
        The given monomial ordering $(kws.ordering) is not implemented for
        $(kws.monoms) monomial representation. Falling back to other monomial
        representations."""
    end
    # if packed representation is requested
    if kws.monoms === :packed
        if is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, ExponentSize}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, ExponentSize}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, ExponentSize}
            end
            @log level = 1 """
            Unable to use $(kws.monoms) monomial representation, too many
            variables ($nvars). Falling back to dense monomial
            representation."""
        else
            @log level = 1 """
            The given monomial ordering $(kws.ordering) is not implemented for
            $(kws.monoms) monomial representation. Falling back to dense
            representation."""
        end
    end
    # in the automatic choice, we always prefer packed representations
    if kws.monoms === :auto
        # TODO: also check that ring.ord is supported
        if is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, ExponentSize}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, ExponentSize}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, ExponentSize}
            end
        end
    end

    ExponentVector{ExponentSize}
end

function get_tight_signed_int_type(x::T) where {T <: Integer}
    if x <= typemax(Int8)
        return Int8
    elseif typemax(Int8) < x <= typemax(Int16)
        return Int16
    elseif typemax(Int16) < x <= typemax(Int32)
        return Int32
    elseif typemax(Int32) < x <= typemax(Int64)
        return Int64
    elseif x <= typemax(Int128)
        return Int128
    else
        @unreachable
        return Int64
    end
end

function get_tight_unsigned_int_type(x::T) where {T <: Integer}
    if x <= typemax(UInt8)
        return UInt8
    elseif typemax(UInt8) < x <= typemax(UInt16)
        return UInt16
    elseif typemax(UInt16) < x <= typemax(UInt32)
        return UInt32
    elseif typemax(UInt32) < x <= typemax(UInt64)
        return UInt64
    elseif x <= typemax(UInt128)
        return UInt128
    else
        @unreachable
        return Int64
    end
end

function select_coefftype(char, npolys, nvars, ordering, kws, hint)
    @log level = -1 """
    Selecting coefficient representation.
    Given hint hint=$hint. Keyword argument arithmetic=$(kws.arithmetic)"""

    if iszero(char)
        return Rational{BigInt}, true
    end
    @assert char > 0

    if char > typemax(UInt64)
        __throw_input_not_supported(
            char,
            "The coefficient field characteristic is too large."
        )
    end

    using_wide_type_for_coeffs = true
    tight_signed_type = get_tight_signed_int_type(char)

    # If the requested arithmetic requires a signed representation
    if kws.arithmetic === :signed
        if typemax(Int32) < char < typemax(UInt32) ||
           typemax(Int64) < char < typemax(UInt64)
            @log level = 1_000 "Cannot use $(kws.arithmetic) arithmetic with characteristic $char"
            @assert false
        elseif !using_wide_type_for_coeffs
            return tight_signed_type, using_wide_type_for_coeffs
        else
            return widen(tight_signed_type), using_wide_type_for_coeffs
        end
    end

    tight_unsigned_type = get_tight_unsigned_int_type(char)
    if !using_wide_type_for_coeffs
        tight_unsigned_type
    else
        widen(tight_unsigned_type)
    end, using_wide_type_for_coeffs
end

###
# Converting to selected polynomial representation

"""
    convert_to_internal(representation, polynomials, kws)

Converts elements of the given array `polynomials` into an internal polynomial
representation specified by the given `representation`.
"""
@timeit function convert_to_internal(
    representation::PolynomialRepresentation,
    polynomials,
    kws::KeywordsHandler;
    dropzeros=true
)
    check_input(polynomials, kws)
    # NOTE: Input polynomials must not be modified.
    @log level = -2 "Converting input polynomials to internal representation.."
    ring = extract_ring(polynomials)
    var_to_index, monoms, coeffs = extract_polys(representation, ring, polynomials)
    @log level = -2 "Done converting input polynomials to internal representation."
    @log level = -6 """
    Polynomials in internal representation:
    Ring: $ring
    Variable to index map: $var_to_index
    Monomials: $monoms
    Coefficients: $coeffs"""
    if dropzeros
        @log level = -2 "Removing zero polynomials"
        remove_zeros_from_input!(ring, monoms, coeffs)
    end
    ring, var_to_index, monoms, coeffs
end

function check_input(polynomials, kws)
    if isempty(polynomials)
        __throw_input_not_supported(polynomials, "Empty input array.")
    end
    _check_input(polynomials, kws)
end

function extract_polys(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    polynomials::Vector{T}
) where {T}
    coeffs = extract_coeffs(representation, ring, polynomials)
    reversed_order, var_to_index, monoms = extract_monoms(representation, ring, polynomials)
    @assert length(coeffs) == length(monoms)
    if reversed_order
        for i in 1:length(coeffs)
            @assert length(coeffs[i]) == length(monoms[i])
            reverse!(coeffs[i])
            reverse!(monoms[i])
        end
    end
    var_to_index, monoms, coeffs
end

function remove_zeros_from_input!(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    @assert length(monoms) == length(coeffs)
    filter!(!iszero_coeffs, coeffs)
    filter!(!iszero_monoms, monoms)
    @assert length(monoms) == length(coeffs)
    iszerobasis = isempty(monoms)
    @log level = -7 "After removing zero polynomials from input:" monoms coeffs
    iszerobasis
end

# Checks that the monomial orderings specified by the given `ring` and
# `params.target_ord` are consistent with the given input monomials `monoms`. In
# case the target ordering differs from the `ring` ordering,  
# sorts the polynomials terms w.r.t. the target ordering.
#
# Also returns the sorting permutations for polynomial terms
function set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    current_ord = ring.ord
    target_ord = params.target_ord
    internal_ord = convert_to_internal_monomial_ordering(var_to_index, target_ord)
    @log level = -2 "Internal ordering:\n$internal_ord"
    ring = PolyRing(ring.nvars, internal_ord, ring.ch)
    if current_ord == target_ord
        # No reordering of terms needed, the terms are already ordered according
        # to the requested monomial ordering
        return ring, Vector{Vector{Int}}()
    end
    @log level = -2 "Reordering input polynomial terms from $(current_ord) to $(target_ord)"
    permutations = sort_input_terms_to_change_ordering!(monoms, coeffs, internal_ord)
    @log level = -7 "Reordered terms:" monoms coeffs
    ring, permutations
end

###
# Converting polynomials from internal representation to the original types

iszero_coeffs(v) = isempty(v)
iszero_monoms(v) = isempty(v)

zero_coeffs(::Type{T}, ring::PolyRing) where {T} = Vector{T}()
zero_monoms(::Type{T}, ring::PolyRing) where {T} = Vector{T}()

"""
    convert_to_output(ring, polynomials, monoms, coeffs, params)

Converts polynomials in internal representation given by arrays `monoms` and
`coeffs` into polynomials in the output format (using `polynomials` as a
reference).
"""
@timeit function convert_to_output(
    ring::PolyRing,
    polynomials,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @assert !isempty(polynomials)
    @log level = -2 "Converting polynomials from internal representation to output format"
    # NOTE: Internal polynomials must not be modified.
    if isempty(monoms)
        @log level = -7 "Output is empty, appending an empty placeholder polynomial"
        push!(monoms, zero_monoms(M, ring))
        push!(coeffs, zero_coeffs(C, ring))
    end
    _convert_to_output(ring, polynomials, monoms, coeffs, params)
end

###
# Utilities for composite coefficients

function unpack_composite_coefficients(
    composite_coeffs::Vector{Vector{CompositeInt{2, T}}}
) where {T <: CoeffZp}
    coeffs_part_1 = Vector{Vector{T}}(undef, length(composite_coeffs))
    coeffs_part_2 = Vector{Vector{T}}(undef, length(composite_coeffs))
    # TODO: Transpose this loop
    @inbounds for i in 1:length(composite_coeffs)
        coeffs_part_1[i] = Vector{T}(undef, length(composite_coeffs[i]))
        coeffs_part_2[i] = Vector{T}(undef, length(composite_coeffs[i]))
        for j in 1:length(composite_coeffs[i])
            a1, a2 = unpack_composite_integer(composite_coeffs[i][j])
            coeffs_part_1[i][j] = a1
            coeffs_part_2[i][j] = a2
        end
    end
    coeffs_part_1, coeffs_part_2
end

function unpack_composite_coefficients(
    composite_coeffs::Vector{Vector{CompositeInt{N, T}}}
) where {N, T <: CoeffZp}
    coeffs_part_i = ntuple(_ -> Vector{Vector{T}}(undef, length(composite_coeffs)), N)
    @inbounds for i in 1:length(composite_coeffs)
        for k in 1:N
            coeffs_part_i[k][i] = Vector{T}(undef, length(composite_coeffs[i]))
        end
        for j in 1:length(composite_coeffs[i])
            ai = unpack_composite_integer(composite_coeffs[i][j])
            for k in 1:N
                coeffs_part_i[k][i][j] = ai[k]
            end
        end
    end
    coeffs_part_i
end
