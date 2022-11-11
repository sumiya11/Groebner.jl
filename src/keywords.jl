# The game with monomial representations:
#  - there are our internal representations (<:RepresentationStyle), 
#  - and user-given hints in the input (<:RepresentationHint).
# For a user-given hint we select the most suitable internal representation.
#
# If the user hint is feasible, we follow it.
# Otherwise, we select some other representation automatically.

# Two monomial representation styles are possible: safe and unsafe
# - A monomal representation is considered safe
#   if a monomial in this representation can safely store a very
#   large exponent.
# - If a monomial cannot store a very large exponent, 
#   then the representation is unsafe.
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
best() = Packed{UInt8}()
# best hint that is also a Safe representation
bestsafe() = NotPacked{UInt64}()
@assert issafe(bestsafe())

# if the given representation is packed, return true if it is good 
# (which just means that there are more than 1 integers in one packed chunk)
isgoodpacked(::Any) = true
isgoodpacked(::Packed{E}) where {E} = true
isgoodpacked(::Packed{UInt64}) = false
isgoodpacked(::Packed{UInt128}) = false

function _not_effective_repr(r)
    @warn "Provided representation hint ($r) looks not effective for input polynomials and was ignored, sorry."
end

default_safe_representation() = default_safe_representation(bestsafe())
function default_safe_representation(hint::T) where {T <: RepresentationHint{E}} where {E}
    @assert issafe(hint)
    Representation{PowerVector{E}}()
end

guess_effective_representation(polynomials) = guess_effective_representation(polynomials, SafeRepresentation())
guess_effective_representation(polynomials, s::RepresentationStyle) = guess_effective_representation(polynomials, s, bestsafe())

# guess effective Safe representation
function guess_effective_representation(polynomials, s::SafeRepresentation, hint::T) where {T<:RepresentationHint{E}} where {E}
    if !issafe(hint) || !isgoodpacked(hint)
        _not_effective_repr(hint)
        hint = bestsafe()
    end
    default_safe_representation(hint)
end

# guess effective Unsafe Notpacked representation
function guess_effective_representation(
        polynomials, 
        s::UnsafeRepresentation, 
        hint::NotPacked{E}) where {E}
    Representation{PowerVector{E}}()
end

# guess effective Unsafe packed representation
function guess_effective_representation(
        polynomials, 
        s::UnsafeRepresentation, 
        hint::Packed{E}) where {E}
    if !isgoodpacked(hint)
        _not_effective_repr(hint)
        hint = best()
    end
    first_impression = peek_at_polynomials(polynomials)
    if first_impression.nvars < div(8, sizeof(E))
        Representation{PackedPair1{UInt64, E}}()
    elseif first_impression.nvars < 2*div(8, sizeof(E))
        Representation{PackedPair2{UInt64, E}}()
    elseif first_impression.nvars < 3*div(8, sizeof(E))
        Representation{PackedPair3{UInt64, E}}()
    else
        Representation{PowerVector{E}}()
    end
end

#------------------------------------------------------------------------------

safe_linear_algebra() = :exact

#------------------------------------------------------------------------------

# Here we choose parameters for groebner basis computation
# based on the specified input keywords

struct GroebnerMetainfo{Rng}
    # if set, then use fglm algorithm for order conversion
    usefglm::Bool
    # output polynomials monomial order
    targetord::Symbol
    # monomial order for computation
    computeord::Symbol

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
    targetord = :lex
    computeord = :lex

    if forsolve
        targetord = :lex
        usefglm = true
        if ordering in (:deglex, :degrevlex)
            computeord = ordering
        else
            computeord = :degrevlex
        end
        computeord = :lex
    else
        if ordering === :input
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
