# Input-output conversions of polynomials.

# Our conventions on some edge cases:
# - Trying to compute a Groebner basis of an empty set is an error.
# - The Groebner basis of [0] is [0].
# - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]

# Polynomials from frontend packages, such as AbstractAlgebra.jl and
# DynamicPolynomials.jl, may not support some of the orderings supported by
# Groebner.jl, e.g., matrix orderings. We compute the basis in the requested
# ordering in Groebner.jl, and then order the terms in the output according to
# some ordering that is supported by the frontend

"""
    PolyRing

Polynomial ring of computation.
"""
mutable struct PolyRing{Ord <: Union{AbstractMonomialOrdering, AbstractInternalOrdering}}
    # Number of variables
    nvars::Int
    # Monomial ordering
    ord::Ord
    # Characteristic of the coefficient ring
    ch::UInt
end

"""
    PolynomialRepresentation

Internal representation of polynomials.
"""
struct PolynomialRepresentation
    monomtype::Type
    coefftype::Type
end

function select_polynomial_representation(polynomials, kws; hint::Symbol=:none)
    frontend, npolys, char, nvars, ordering = peek_at_polynomials(polynomials)
    monomtype = select_monomtype(char, npolys, nvars, kws, hint)
    coefftype = select_coefftype(char, npolys, nvars, kws, hint)
    basering = iszero(char) ? :qq : :zp
    @log level = -1 "Frontend: $frontend"
    @log level = -1 """
    Input: $npolys polynomials over $basering in $nvars variables
    Ordering: $ordering"""
    @log level = -1 """
    Internal representation: 
    monomials are $monomtype
    coefficients are $coefftype"""
    PolynomialRepresentation(monomtype, coefftype)
end

function select_monomtype(char, npolys, nvars, kws, hint)
    # TODO: also dispatch on the type of monomial ordering    
    @log level = -1 "Selecting monomial representation.\nGiven hint hint=$hint. Keyword argument monoms=$(kws.monoms)"
    if hint === :large_exponents
        @log level = -1 "As hint=$hint was provided, using 64 bits per single exponent"
        # use 64 bits if large exponents detected
        @assert is_supported_ordering(ExponentVector{UInt64}, kws.ordering)
        return ExponentVector{UInt64}
    end
    E = UInt8
    variables_per_word = div(sizeof(UInt), sizeof(E))
    # if dense representation is requested
    if kws.monoms === :dense
        @assert is_supported_ordering(ExponentVector{E}, kws.ordering)
        return ExponentVector{E}
    end
    # if sparse representation is requested
    if kws.monoms === :sparse
        if is_supported_ordering(SparseExponentVector{E, nvars, Int}, kws.ordering)
            return SparseExponentVector{E, nvars, Int}
        end
        @log level = 1 """
        The given monomial ordering $(kws.ordering) is not implemented for $(kws.monoms) monomials, sorry.
        Falling back to dense representation."""
    end
    # if packed representation is requested
    if kws.monoms === :packed
        if is_supported_ordering(PackedTuple1{UInt64, E}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, E}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, E}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, E}
            end
            @log level = 1 """
            Unable to use $(kws.monoms) monomials, too many variables ($nvars), sorry.
            Falling back to dense monomial representation."""
        else
            @log level = 1 """
            The given monomial ordering $(kws.ordering) is not implemented for $(kws.monoms) monomials, sorry.
            Falling back to dense representation."""
        end
    end
    if kws.monoms === :auto
        if is_supported_ordering(PackedTuple1{UInt64, E}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, E}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, E}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, E}
            end
        end
    end
    ExponentVector{E}
end

function select_coefftype(char, npolys, nvars, kws, hint)
    if !iszero(char)
        @assert char < typemax(UInt64)
        if char > 2^32
            UInt128
        else
            UInt64
        end
    else
        Rational{BigInt}
    end
end

"""
    convert_to_internal(polynomials, kws, representation)

Converts elements of the given array `polynomials` into an internal polynomial
representation specified by `representation`.

Returns a tuple (`ring`, `monoms`, `coeffs`).
"""
function convert_to_internal(
    representation::PolynomialRepresentation,
    polynomials,
    kws::KeywordsHandler;
    dropzeros=true
)
    check_input(polynomials)
    # NOTE: internally, 0 polynomial is represented with an empty vector of
    # monomials and an empty vector of coefficients.
    # NOTE: Input polynomials must not be modified.
    @log level = -2 "Converting input polynomials to internal representation.."
    ring = extract_ring(polynomials)
    # @assert is_representation_suitable(representation, ring)
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

function check_input(polynomials)
    !isempty(polynomials) && "Empty input polynomials"
    # TODO: check that there are no parameters
    # TODO: check that the ground field is Z_p or QQ
    nothing
end

function extract_polys(representation, ring::PolyRing, polynomials::Vector{T}) where {T}
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

"""
    convert_to_output(ring, polynomials, monoms, coeffs, params)

Converts internal polynomials given by arrays `monoms` and `coeffs` into
polynomials in the output format.
"""
function convert_to_output(
    ring,
    polynomials,
    monoms::Vector{M},
    coeffs::Vector{C},
    params
) where {M, C}
    @assert !isempty(polynomials)
    @log level = -2 "Converting polynomials from internal representation to output format"
    # NOTE: Internal polynomials must not be modified.
    if isempty(monoms)
        push!(monoms, Vector{M}())
        push!(coeffs, Vector{C}())
    end
    origring = parent(first(polynomials))
    convert_to_output(origring, monoms, coeffs, params)
end

# checks that the coefficient `c` can be represented exactly in type `T`.
checkexact(c, T::Type{BigInt}) = true
checkexact(c, T::Type{Rational{U}}) where {U} =
    checkexact(numerator(c), U) && checkexact(denominator(c), U)
checkexact(c, T) = typemin(T) <= c <= typemax(T)

@noinline function __throw_inexact_coeff_conversion(c, T)
    throw(
        DomainError(
            c,
            """
            Coefficient $c in the output basis cannot be converted exactly to $T. 
            Using big arithmetic in the input should fix this."""
        )
    )
end

function check_and_convert_coeffs(coeffs_zz, T)
    cfs = Vector{T}(undef, length(coeffs_zz))
    for i in 1:length(coeffs_zz)
        !checkexact(coeffs_zz[i], T) && __throw_inexact_coeff_conversion(coeffs_zz[i], T)
        cfs[i] = coeffs_zz[i]
    end
    cfs
end

function convert_coeffs_to_output(
    coeffs::Vector{Q},
    ::Type{T}
) where {Q <: CoeffQQ, T <: Rational}
    check_and_convert_coeffs(coeffs, T)
end

function convert_coeffs_to_output(
    coeffs::Vector{Q},
    ::Type{T}
) where {Q <: CoeffQQ, T <: Integer}
    coeffs_zz = clear_denominators(coeffs)
    check_and_convert_coeffs(coeffs_zz, T)
end

iszero_coeffs(v) = isempty(v)
iszero_monoms(v) = isempty(v)

zero_coeffs_ff(::Type{T}, ring::PolyRing) where {T} = T[]
zero_coeffs_qq(::Type{T}, ring::PolyRing) where {T} = T[]

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
function set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    current_ord = ring.ord
    target_ord = params.target_ord
    internal_ord = convert_to_internal_monomial_ordering(var_to_index, target_ord)
    @log level = -2 "Internal ordering:\n$internal_ord"
    ring = PolyRing(ring.nvars, internal_ord, ring.ch)
    if current_ord == target_ord
        # No reordering of terms needed, they are already ordered according to
        # the requested ordering
        return ring
    end
    @log level = -2 "Reordering input polynomial terms from $(current_ord) to $(target_ord)"
    sort_input_to_change_ordering!(monoms, coeffs, internal_ord)
    ring
end
