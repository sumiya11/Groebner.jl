# The game with monomial representations:
#  - there are our internal representations (<:RepresentationStyle), 
#  - and user-given hints in the input (<:RepresentationHint).
# For a user-given hint we select the most suitable internal representation.
#
# If the user hint is feasible, we follow it.
# Otherwise, we select some other representation automatically.

_not_effective_repr(r) = @warn "Provided monom representation hint ($r) was ignored, sorry."

# Two monomial representation styles are possible: safe and unsafe
# - A monomal representation is considered safe
#   if a monomial in this representation can safely store a very large exponent.
# - If a monomial cannot store a very large exponent, 
#   then the representation is considered unsafe.
abstract type RepresentationStyle end
struct SafeRepresentation <: RepresentationStyle end
struct UnsafeRepresentation <: RepresentationStyle end

struct Representation{M<:Monom} <: RepresentationStyle end

# The number of variables in monomial representation max.
capacity(::Representation{M}) where {M<:Monom} = capacity(M)

# The user can give a hint of the desired monomal representation,
# Possible hints are:
# - Packed{E<:Unsigned} - packed representation with sizeof(E)*8 bits per exponent 
# - NotPacked{E<:Unsigned} - not packed representation with sizeof(E)*8 bits per exponent
# - Sparse{E<:Unsigned} (currently not used)
abstract type RepresentationHint{E} end
struct NotPacked{E<:Unsigned} <: RepresentationHint{E} end
struct Packed{E<:Unsigned} <: RepresentationHint{E} end
struct Sparse{E<:Unsigned} <: RepresentationHint{E} end

# representation is safe if it can store exponents up to at least 2^32 - 1
issafe(::T) where {T <: RepresentationHint{E}} where {E} = sizeof(E) >= 4

# best hint (should be default)
best_monom_representation() = Packed{UInt8}()
# best hint that is also a Safe representation
bestsafe() = NotPacked{UInt64}()
@assert issafe(bestsafe())

# if the given representation is packed, return true if it is a good packed representation
# (which just means that there are more than 1 integers in one packed chunk)
isgoodpacked(::Any) = true
isgoodpacked(::Packed{E}) where {E} = true
isgoodpacked(::Packed{UInt64}) = false
isgoodpacked(::Packed{UInt128}) = false

default_safe_representation() = default_safe_representation(bestsafe())
function default_safe_representation(hint::T) where {T <: RepresentationHint{E}} where {E}
    @assert issafe(hint)
    Representation{PowerVector{E}}()
end

guess_effective_representation(polynomials) = guess_effective_representation(polynomials, SafeRepresentation())
guess_effective_representation(polynomials, s::RepresentationStyle) = guess_effective_representation(polynomials, s, bestsafe())

# guess effective Safe representation
function guess_effective_representation(
        polynomials, s::SafeRepresentation, ordering, hint::T
        ) where {T<:RepresentationHint{E}} where {E}
    if !issafe(hint) || !isgoodpacked(hint)
        _not_effective_repr(hint)
        hint = bestsafe()
    end
    default_safe_representation(hint)
end

# guess effective Unsafe Not packed representation
function guess_effective_representation(
        polynomials, 
        s::UnsafeRepresentation, 
        ordering, 
        hint::NotPacked{E}) where {E}
    @assert is_supported_ordering(PowerVector{E}, ordering)
    Representation{PowerVector{E}}()
end

# guess effective Unsafe packed representation
function guess_effective_representation(
        polynomials, 
        s::UnsafeRepresentation, 
        ordering, 
        hint::Packed{E}) where {E}
    if !isgoodpacked(hint)
        _not_effective_repr(hint)
        hint = best_monom_representation()
    end
    first_impression = peek_at_polynomials(polynomials)
    elper8bytes = div(8, sizeof(E))
    # if we want a non-packed representation
    if first_impression.nvars > 3*elper8bytes
        @assert is_supported_ordering(PowerVector{E}, ordering)
        return Representation{PowerVector{E}}()
    end
    # if we want a packed representation
    if is_supported_ordering(AbstractPackedPair, ordering)
        if first_impression.nvars < elper8bytes
            return Representation{PackedPair1{UInt64, E}}()
        elseif first_impression.nvars < 2*elper8bytes
            return Representation{PackedPair2{UInt64, E}}()
        elseif first_impression.nvars < 3*elper8bytes
            return Representation{PackedPair3{UInt64, E}}()
        end
    end
    # return the default safe representation
    return default_safe_representation()
end

function peek_at_polynomials(polynomials::Vector{T}) where {T}
    (nvars=2^32,)
end

function peek_at_polynomials(polynomials::Vector{T}) where {T<:AbstractAlgebra.Generic.MPolyElem}
    isempty(polynomials) && return (nvars=2^32,)
    (
        nvars=AbstractAlgebra.nvars(parent(polynomials[1])),
    )
end

#------------------------------------------------------------------------------

# Default safe option of linear algebra 
safe_linear_algebra() = :exact

#------------------------------------------------------------------------------

# Here we choose parameters for groebner basis computation
# based on the specified input keywords

struct GroebnerMetainfo{Rng, Ord1, Ord2}
    # if set, then use fglm algorithm for order conversion
    usefglm::Bool
    # output polynomials monomial order
    targetord::Ord1
    # monomial order for computation
    computeord::Ord2

    # correctness checks levels
    heuristiccheck::Bool
    randomizedcheck::Bool
    guaranteedcheck::Bool

    # linear algebra backend to be used
    # Currently available are
    #   :exact for exact linear algebra,
    #   :prob for probabilistic linear algebra
    linalg::Symbol

    # ground field of computation.
    # Current options are
    #   :qq for rationals,
    #   :ff for integers modulo prime
    ground::Symbol

    # random number generator
    rng::Rng
end

function set_metaparameters(ring, ordering, certify, forsolve, linalg, rng)
    usefglm = false
    targetord = Lex()
    computeord = Lex()

    if forsolve
        targetord = Lex()
        usefglm = true
        if ordering in (DegLex(), DegRevLex())
            computeord = ordering
        else
            computeord = DegRevLex()
        end
        computeord = Lex()
    else
        if ordering === InputOrdering()
            ordering = ring.ord
        end
        targetord = ordering
        usefglm = false
        computeord = targetord
    end

    heuristiccheck = true
    randomizedcheck = true
    if certify
        guaranteedcheck = true
    else
        guaranteedcheck = false
    end

    ground = :qq
    if ring.ch == 0
        ground = :qq
    else
        ground = :ff
    end

    @info "Computing in $computeord order, result is in $targetord order"
    @info "Using fglm: $usefglm"

    GroebnerMetainfo(usefglm, targetord, computeord,
                        heuristiccheck, randomizedcheck, guaranteedcheck,
                        linalg, ground, rng)
end
