# Input-output conversions of polynomials.
# Conversions between input polynomials and internal polynomial representations

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

    Information about polynomial ring.
"""
mutable struct PolyRing{Ord <: AbstractMonomialOrdering}
    # Number of variables
    nvars::Int
    # Monomial ordering
    ord::Ord
    # Characteristic of the coefficient ring
    ch::Int
end

"""
    PolynomialRepresentation


"""
struct PolynomialRepresentation
    monomtype::Type
    coefftype::Type
end

function determine_monomtype(basering, npolys, nvars, kws)
    E = UInt8
    elper8bytes = div(8, sizeof(E))
    # if we want a non-packed representation
    if nvars > 3 * elper8bytes
        @assert is_supported_ordering(PowerVector{E}, kws.ordering)
        return PowerVector{E}
    end
    # if we want a packed representation
    if is_supported_ordering(AbstractPackedPair, kws.ordering)
        if nvars < elper8bytes
            return PackedPair1{UInt64, E}
        elseif nvars < 2 * elper8bytes
            return PackedPair2{UInt64, E}
        elseif nvars < 3 * elper8bytes
            return PackedPair3{UInt64, E}
        end
    end
    PowerVector{E}
end

function determine_coefftype(basering, npolys, nvars, kws)
    basering == :zp ? UInt64 : Rational{BigInt}
end

function select_polynomial_representation(polynomials, kws; hint=nothing)
    frontend, npolys, basering, nvars, ordering = peek_at_polynomials(polynomials)
    monomtype = determine_monomtype(basering, npolys, nvars, kws)
    coefftype = determine_coefftype(basering, npolys, nvars, kws)
    @log level = 1 "Frontend: $frontend"
    @log level = 1 "Input: $npolys polynomials over $basering in $nvars variables, $ordering ordering"
    @log level = 1 "Internal representation: monomials are $monomtype, coefficients are $coefftype"
    PolynomialRepresentation(monomtype, coefftype)
end

"""
    peek_at_polynomials(polynomials)

Takes a peek at input polynomials and returns a tuple 
(frontend, npolys, basering, nvars, ordering)
"""
function peek_at_polynomials end

"""
    convert_to_internal(polynomials, kws, representation)

Converts elements of an array `polynomials` into an internal polynomial
representation specified by `representation`.

Returns a tuple (`allzeros`, `ring`, `monoms`, `coeffs`).
"""
function convert_to_internal(
    representation::PolynomialRepresentation,
    polynomials,
    kws::KeywordsHandler
)
    check_input(polynomials)
    # NOTE: internally, 0 polynomial is represented with an empty vector of
    # monomials and an empty vector of coefficients.
    # NOTE: Input polynomials must not be modified.
    @log level = 3 "Converting input polynomials to internal representation.."
    ring = extract_ring(polynomials)
    # @assert is_representation_suitable(representation, ring)
    monoms, coeffs = extract_polys(representation, ring, polynomials)
    # allzeros = remove_zeros_from_input!(ring, monoms, coeffs)
    @log level = 3 "Done converting input polynomials to internal representation."
    @log level = -100 """
    Polynomials in internal representation:
    Ring: $ring
    Monomials: $monoms
    Coefficients: $coeffs"""
    ring, monoms, coeffs
end

function check_input(polynomials)
    !isempty(polynomials) && "Empty input"
    # TODO: check that there are no parameters
    # TODO: check that the ground field is Z_p or QQ
    nothing
end

function extract_polys(representation, ring::PolyRing, polynomials::Vector{T}) where {T}
    coeffs = extract_coeffs(representation, ring, polynomials)
    monoms = extract_monoms(representation, ring, polynomials)
    monoms, coeffs
end

"""
    convert_to_output(ring, polynomials, monoms, coeffs, params)

Converts internal polynomials given by arrays `monoms` and `coeffs` into
polynomials in the output format.
"""
function convert_to_output(ring, polynomials, monoms, coeffs, params)
    # NOTE: Internal polynomials must not be modified.
    # TODO: throw warning if the output format is strange
    convert_to_output(parent(first(polynomials)), monoms, coeffs, params)
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
            """Coefficient $c in the output basis cannot be converted exactly to $T. 
            Using big arithmetic in the input should fix this.
            """
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
    coeffs_zz = scale_denominators(coeffs)
    check_and_convert_coeffs(coeffs_zz, T)
end

iszero_coeffs(v) = isempty(v)
iszero_monoms(v) = isempty(v)

zero_coeffs_ff(ring::PolyRing{Ch}) where {Ch} = Ch[]
zero_coeffs_qq(ring::PolyRing{Ch}) where {Ch} = Rational{BigInt}[]

function remove_zeros_from_input!(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    @assert length(monoms) == length(coeffs)
    filter!(iszero_coeffs, coeffs)
    filter!(iszero_monoms, monoms)
    @assert length(monoms) == length(coeffs)
    iszerobasis = isempty(monoms)
    if iszerobasis
        push!(monoms, Vector{M}())
        push!(coeffs, Vector{T}())
    end
    iszerobasis
end

# Check the consistency of the given monomial ordering
# with respect to the monomial implementation and the length of exponent vector

# Should this be moved to src/monoms ?
@noinline function __throw_monomial_ordering_inconsistent(e, o)
    throw(
        DomainError(
            o,
            """The given monomial ordering is inconsistent with the input.
            Monomial: $e
            Ordering: $o
            Possible cause is that the number of variables in the ordering is set incorrectly."""
        )
    )
end

check_ordering(e::M, o::Lex) where {M <: Monom} = true
check_ordering(e::M, o::DegLex) where {M <: Monom} = true
check_ordering(e::M, o::DegRevLex) where {M <: Monom} = true
function check_ordering(e::M, o::Union{Lex, DegLex, DegRevLex}, lo, hi) where {M <: Monom}
    if lo <= hi
        true
    else
        false
    end
end

function check_ordering(
    e::M,
    wo::WeightedOrdering{O}
) where {M <: Monom, O <: AbstractMonomialOrdering}
    false
end
function check_ordering(
    e::PowerVector{T},
    wo::WeightedOrdering{O},
    lo::Int,
    hi::Int
) where {T, O <: AbstractMonomialOrdering}
    check_ordering(e, wo.ord, lo, hi)
    if hi - lo + 1 != length(wo.weights)
        return false
    end
    true
end
function check_ordering(
    e::PowerVector{T},
    wo::WeightedOrdering{O}
) where {T, O <: AbstractMonomialOrdering}
    check_ordering(e, wo, 2, length(e))
end

function check_ordering(
    e::M,
    bo::BlockOrdering{R1, R2, O1, O2}
) where {M <: Monom, R1, R2, O1 <: AbstractMonomialOrdering, O2 <: AbstractMonomialOrdering}
    false
end
function check_ordering(
    e::PowerVector{T},
    bo::BlockOrdering{R1, R2, O1, O2},
    lo::Int,
    hi::Int
) where {T, R1, R2, O1 <: AbstractMonomialOrdering, O2 <: AbstractMonomialOrdering}
    r1 = (first(bo.r1) + 1):(last(bo.r1) + 1)
    r2 = (first(bo.r2) + 1):(last(bo.r2) + 1)
    if first(r1) != lo || last(r2) != hi
        return false
    end
    check_ordering(e, bo.ord1, first(r1), last(r1))
    check_ordering(e, bo.ord2, first(r2), last(r2))
    true
end
function check_ordering(
    e::PowerVector{T},
    bo::BlockOrdering{R1, R2, O1, O2}
) where {T, R1, R2, O1 <: AbstractMonomialOrdering, O2 <: AbstractMonomialOrdering}
    check_ordering(e, bo, 2, length(e))
end

function check_ordering(e::M, mo::MatrixOrdering) where {M <: Monom}
    false
end
function check_ordering(e::PowerVector{T}, mo::MatrixOrdering) where {T}
    check_ordering(e, mo, 2, length(e))
end
function check_ordering(e::PowerVector{T}, mo::MatrixOrdering, lo::Int, hi::Int) where {T}
    rows = mo.rows
    n = hi - lo + 1
    for i in 1:length(rows)
        if length(rows[i]) != n
            return false
        end
    end
    true
end

# Checks that the monomial orderings specified by the given `ring` and `target_ord` 
# are consistent with the given input monomials `monoms`.
# In case the target ordering differs from the `ring` ordering,  
# sorts the polynomials terms w.r.t. the target ordering.
function change_ordering_if_needed!(
    ring,
    monoms,
    coeffs,
    params
)
    @assert !isempty(monoms) && !isempty(monoms[1])
    !check_ordering(monoms[1][1], params.target_ord) &&
        __throw_monomial_ordering_inconsistent(monoms[1][1], params.target_ord)
    ring.ord == params.target_ord && return nothing
    @log level = 3 "Reordering polynomial terms from $(ring.ord) to $target_ord"
    ring.ord = params.target_ord
    sort_input_to_change_ordering!(monoms, coeffs, params.target_ord)
    nothing
end
