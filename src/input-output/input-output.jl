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
            Using big arithmetic in the input should fix this."""
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

"""
    select_polynomial_representation(polynomials, keywords; hint=:none)

Given an array of input polynomials tries to select a suitable representation
for coefficients and exponent vectors.

Additionally, `hint` can be specified to one of the following:

- `:large_exponents`: use at least 32 bits per exponent.
"""
function select_polynomial_representation(
    polynomials,
    kws::KeywordsHandler;
    hint::Symbol=:none
)
    if !(hint in (:none, :large_exponents))
        @log level = 1000 "The given hint=$hint was discarded"
    end
    frontend, npolys, char, nvars, ordering = peek_at_polynomials(polynomials)
    monomtype = select_monomtype(char, npolys, nvars, ordering, kws, hint)
    coefftype = select_coefftype(char, npolys, nvars, ordering, kws, hint)
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

function select_monomtype(char, npolys, nvars, ordering, kws, hint)
    @log level = -1 "Selecting monomial representation.\nGiven hint hint=$hint. Keyword argument monoms=$(kws.monoms)"
    if hint === :large_exponents
        @log level = -1 "As hint=$hint was provided, using 64 bits per single exponent"
        # use 64 bits if large exponents detected
        desired_monom_type = ExponentVector{UInt64}
        @assert is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
    end
    # TODO: explain this condition
    if kws.homogenize === :yes || (
        kws.homogenize === :auto && (
            kws.ordering isa Lex ||
            kws.ordering isa ProductOrdering ||
            (ordering == :lex && kws.ordering isa InputOrdering)
        )
    )
        @log level = -1 "As homogenize=:yes/:auto was provided, representing monomials as exponent vectors"
        desired_monom_type = ExponentVector{UInt32}
        @assert is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
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
        The given monomial ordering $(kws.ordering) is not implemented for $(kws.monoms) monomials.
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
            Unable to use $(kws.monoms) monomials, too many variables ($nvars).
            Falling back to dense monomial representation."""
        else
            @log level = 1 """
            The given monomial ordering $(kws.ordering) is not implemented for $(kws.monoms) monomials.
            Falling back to dense representation."""
        end
    end
    if kws.monoms === :auto
        # TODO: also check that ring.ord is supported
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

function select_coefftype(char, npolys, nvars, ordering, kws, hint)
    if !iszero(char)
        if char >= typemax(UInt64)
            __throw_input_not_supported(char, "The coefficient field order is too large.")
        end
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
    convert_to_internal(representation, polynomials, kws)

Converts elements of the given array `polynomials` into an internal polynomial
representation specified by `representation`.

Returns a tuple (`ring`, `var_to_index`, `monoms`, `coeffs`).
"""
function convert_to_internal(
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
    # TODO: check that there are no parameters
    # TODO: check that the ground field is Z_p or QQ
    true
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
function convert_to_output(
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
