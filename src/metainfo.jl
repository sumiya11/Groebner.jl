
abstract type RepresentationStyle end
struct SafeRepresentation <: RepresentationStyle end
struct UnsafeRepresentation <: RepresentationStyle end

struct Representation{M<:Monom} <: RepresentationStyle end

capacity(::Representation{M}) where {M<:Monom} = capacity(M)

abstract type RepresentationHint{E} end
struct NotPacked{E<:Unsigned} <: RepresentationHint{E} end
struct Packed{E<:Unsigned} <: RepresentationHint{E} end
struct Sparse{E<:Unsigned} <: RepresentationHint{E} end

best() = Packed{UInt8}()
bestsafe() = NotPacked{UInt64}()
issafe(::T) where {T <: RepresentationHint{E}} where {E} = sizeof(E) >= 4
@assert issafe(bestsafe())

isgoodpacked(::Any) = true
isgoodpacked(::Packed{E}) where {E} = true
isgoodpacked(::Packed{UInt64}) = false

function _not_effective_repr(r)
    @warn "Provided monomial representation ($r) looks not effective for input polynomials and was ignored, sorry."
end

default_safe_representation() = default_safe_representation(bestsafe())
function default_safe_representation(hint::T) where {T <: RepresentationHint{E}} where {E}
    @assert issafe(hint)
    Representation{PowerVector{E}}()
end

guess_effective_representation(polynomials) = guess_effective_representation(polynomials, SafeRepresentation())
guess_effective_representation(polynomials, s::RepresentationStyle) = guess_effective_representation(polynomials, s, bestsafe())

function guess_effective_representation(polynomials, s::SafeRepresentation, hint::T) where {T<:RepresentationHint{E}} where {E}
    if !issafe(hint) || !isgoodpacked(hint)
        _not_effective_repr(hint)
        hint = bestsafe()
    end
    default_safe_representation(hint)
end

function guess_effective_representation(
        polynomials, 
        s::UnsafeRepresentation, 
        hint::NotPacked{E}) where {E}
    Representation{PowerVector{E}}()
end

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
    #   :exact for exact linear algebra
    #   :prob for probabilistic linear algebra
    linalg::Symbol
    
    ground::Symbol

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
        # TODO: lex fglm
        computeord = :lex
    else
        if ordering == :input
            ordering = ring.ord
        end
        targetord = ordering
        usefglm = false
        computeord = targetord
    end

    heuristiccheck = true
    # heuristiccheck = false
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
