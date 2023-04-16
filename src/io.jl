# Input-output conversions of polynomials.
# Currently, conversions work with polynomials from
#  - AbstractAlgebra.jl
#  - DynamicPolynomials.jl
#  - Nemo.jl
#  - Singular.jl (TODO: add tests)

#=
    Our conventions:
    - Trying to compute a Groebner basis of an empty set is an error.
    - The Groebner basis of [0] is [0].
    - The Groebner basis of [0,..., 0] is [0].
    - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]
=#

#=
    A note on how we represent polynomials.

    First, coefficients, exponents and polynomial ring information
    are extracted from input polynomials with `convert_to_internal`.

    Inside the algorithm all monomials are hashed without collisions,
    so that an integer represents a single monomial.
    A single monomial could be stored in several different ways, 
    according to the prescibed monomial implementation.

    A polynomial is represented with a dynamic coefficients vector
    together with a dynamic vector of hashtable indices of monomials.

    After the basis is computed, hash table indices are converted back to 
    monomials, and the `convert_to_output` is called 
    to convert internal structures to the original polynomial type.
=#

#=
    A note about monomial orderings.

    Groebner.jl supports several monomial orderings 
    (all of which are subtypes of AbstractMonomialOrdering, 
    see src/monoms/orderings.jl for details)

    Polynomials from AbstractAlgebra.jl, DynamicPolynomials.jl, 
    and some other packages, do not support some of the orderings supported by Groebner.jl.

    We compute the basis in the requested ordering in Groebner.jl,
    and then output the polynomials in the ordering of the input.

    For example, say that the input polynomials are from AbstractAlgebra.jl 
    and are in the Lex ordering. The requested basis is in some Weighted ordering.
    Then, the output should be a correct basis in the Weighted ordering,
    though the terms in the output are ordered w.r.t. Lex ordering.
=#

# For compatibility with polynomials from AbstractAlgebra.jl and Nemo.jl
const _AA_supported_orderings_symbols = (:lex, :deglex, :degrevlex)
# The type of exponent vector entries used internally in AbstractAlgebra.jl
const _AA_exponenttype = UInt64

"""
    Contains information about some polynomial ring.
"""
mutable struct PolyRing{Char<:CoeffFF, Ord<:AbstractMonomialOrdering}
    # number of variables
    nvars::Int
    # monomial ordering
    ord::Ord
    # characteristic of the coefficient field
    ch::Char
    # information about the original ring of input. Options are:
    #    :AbstractAlgebra for AbstractAlgebra,
    #    :MultivariatePolynomials for subtypes of MultivariatePolynomials.jl, 
    #       e.g, for DynamicPolynomials.jl,
    #    :hasparent for polynomials constructed with parent ring, e.g., Nemo
    #    :undefined
    origring::Symbol
end

"""
    Converts input polynomials to internal representation used by the algorithm.
    Extracts polynomial ring information, and polynomial exponents and coefficients.

    This is the most general implementation.
    Works for polynomials that implement AbstractAlgebra.Generic.MPoly intefrace:
        . `AbstractAlgebra`
        . `Nemo`
        . `Singular`

    Currently, there are also specializations for polynomials from:
        . `MultivariatePolynomials`

"""
function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::AbstractMonomialOrdering) where {T}
    isempty(orig_polys) && throw(DomainError(orig_polys, "Empty input."))

    if hasmethod(AbstractAlgebra.parent, Tuple{typeof(first(orig_polys))})
        convert_to_internal(representation, orig_polys, ordering, Val(:hasparent))
    else
        convert_to_internal(representation, orig_polys, ordering, Val(:undefined))
    end
end

function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::AbstractMonomialOrdering,
        ::Val{:undefined}) where {T}
    error("Sorry, we don't work with this type of polynomials ($T) yet. Feel free to open an issue")
end

"""
    `hasparent` convention specialization
"""
function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::AbstractMonomialOrdering,
        ::Val{:hasparent}) where {T}
    R = parent(first(orig_polys))
    ring = extract_ring(R)
    check_domain(representation, ring)
    exps, cfs = extract_polys(representation, ring, orig_polys, ring.ord)
    ring, exps, cfs
end

"""
    `MultivariatePolynomials.AbstractPolynomial` conversion specialization
"""
function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::AbstractMonomialOrdering,
        ::Val{:undefined}) where {T<:AbstractPolynomialLike{U}} where {U}
    ring = extract_ring(orig_polys)
    check_domain(representation, ring)
    exps, cfs = extract_polys(representation, ring, orig_polys, ring.ord)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

# Determines the monomial ordering of the output,
# given the original ordering `origord` and the targer ordering `targetord`
ordering_typed2sym(origord, targetord::Lex) = :lex
ordering_typed2sym(origord, targetord::DegLex) = :deglex
ordering_typed2sym(origord, targetord::DegRevLex) = :degrevlex
ordering_typed2sym(origord) = origord
ordering_typed2sym(origord, targetord::AbstractMonomialOrdering) = origord

function ordering_sym2typed(ord::Symbol)
    ord in _AA_supported_orderings_symbols || throw(DomainError(ord, "Not a supported ordering."))
    if ord === :lex
        Lex()
    elseif ord === :deglex
        DegLex()
    elseif ord === :degrevlex
        DegRevLex()
    end
end

function determinechartype(ch)
    Ch = UInt128 
    if ch == 0
        Ch = UInt64
    elseif ch < 2^4
        Ch = UInt8
    elseif ch < 2^8
        Ch = UInt16
    elseif ch < 2^16
        Ch = UInt32
    elseif ch < 2^32
        Ch = UInt64
    end
    Ch
end

function check_ground_domain(nv, ord, ch)
    # @assert ord in _supported_orderings
    @assert 0 <= ch < typemax(UInt64)
end

function check_domain(representation, ring)
    @assert capacity(representation) >= ring.nvars
end

function extract_ring(R::T) where {T}
    @assert hasmethod(AbstractAlgebra.nvars, Tuple{T})
    @assert hasmethod(AbstractAlgebra.characteristic, Tuple{T})

    nv     = AbstractAlgebra.nvars(R)
    # deglex is the default ordering on univariate polynomials
    ord    = hasmethod(AbstractAlgebra.ordering, Tuple{T}) ? AbstractAlgebra.ordering(R) : :deglex
    # type unstable:
    ordT   = ordering_sym2typed(ord)
    ch     = AbstractAlgebra.characteristic(R)

    check_ground_domain(nv, ordT, ch)

    Char = determinechartype(ch)

    PolyRing{Char, typeof(ordT)}(nv, ordT, Char(BigInt(ch)), :AbstractAlgebra)
end

function extract_ring(orig_polys::Vector{<:AbstractPolynomialLike{T}}) where {T}
    f = first(orig_polys)

    nv = Groebner.MultivariatePolynomials.nvariables(orig_polys)
    ordT   = DegLex()
    ch     = 0

    check_ground_domain(nv, ordT, ch)

    Char = determinechartype(ch)

    PolyRing{Char, typeof(ordT)}(nv, ordT, Char(ch), :MultivariatePolynomials)
end

#------------------------------------------------------------------------------

function extract_polys(
        representation,
        ring::PolyRing,
        orig_polys::Vector{T},
        ord::AbstractMonomialOrdering) where {T}
    cfs = extract_coeffs(ring, orig_polys)
    exps = extract_exponents(representation, ring, orig_polys)
    exps, cfs
end

#------------------------------------------------------------------------------

# Our convention is that an empty vector of coefficients
# represents zero polynomial

function iszero_coeffvector(v)
    isempty(v)
end

function iszero_monomvector(v)
    isempty(v)
end

function zero_coeffvector_ff(ring::PolyRing{Ch}) where {Ch}
    Ch[]
end

function zero_coeffvector_qq(ring::PolyRing{Ch}) where {Ch}
    Rational{BigInt}[]
end

function extract_coeffs(ring::PolyRing{Ch}, orig_polys::Vector{T}) where {Ch, T}
    if ring.ch > 0
        extract_coeffs_ff(ring, orig_polys)
    else
        extract_coeffs_qq(ring, orig_polys)
    end
end

# specialization for univariate polynomials
function extract_coeffs_ff(ring::PolyRing{Ch}, poly::Union{AbstractAlgebra.Generic.Poly, AbstractAlgebra.PolyElem}) where {Ch}
    iszero(poly) && (return zero_coeffvector_ff(ring))
    reverse(map(Ch ∘ AbstractAlgebra.data, filter(!iszero, collect(AbstractAlgebra.coefficients(poly)))))
end

# specialization for univariate polynomials
function extract_coeffs_qq(ring::PolyRing, poly::Union{AbstractAlgebra.Generic.Poly, AbstractAlgebra.PolyElem})
    iszero(poly) && (return zero_coeffvector_qq(ring))
    reverse(map(Rational, filter(!iszero, collect(AbstractAlgebra.coefficients(poly)))))
end

# specialization for multivariate polynomials
function extract_coeffs_ff(ring::PolyRing{Ch}, poly) where {Ch}
    iszero(poly) && (return zero_coeffvector_ff(ring))
    map(Ch ∘ AbstractAlgebra.data, AbstractAlgebra.coefficients(poly))
end

# specialization for multivariate polynomials
function extract_coeffs_qq(ring::PolyRing, poly)
    iszero(poly) && (return zero_coeffvector_qq(ring))
    map(Rational, AbstractAlgebra.coefficients(poly))
end

function extract_coeffs_ff(ring::PolyRing{Ch}, 
                    orig_polys::Vector{T}) where {Ch, T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Ch}}(undef, npolys)
    @inbounds for i in 1:npolys
        coeffs[i] = extract_coeffs_ff(ring, orig_polys[i])
    end
    coeffs
end

function extract_coeffs_qq(ring::PolyRing, orig_polys::Vector{T}) where {T}
    map(poly -> extract_coeffs_qq(ring, poly), orig_polys)
end

function extract_coeffs_qq(ring::PolyRing, poly::T) where {T<:AbstractPolynomialLike{U}} where {U}
    iszero(poly) && (return zero_coeffvector_qq(ring))
    map(Rational, MultivariatePolynomials.coefficients(poly))
end

function extract_coeffs_qq(
            ring::PolyRing,
            orig_polys::Vector{T}) where {T<:AbstractPolynomialLike{U}} where {U}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = extract_coeffs_qq(ring, poly)
    end
    coeffs
end

#------------------------------------------------------------------------------

function extract_exponents(representation::Representation{M}, ring::PolyRing, poly) where {M}
    exps = Vector{M}(undef, length(poly))
    @inbounds for j in 1:length(poly)
        exps[j] = make_ev(M, AbstractAlgebra.exponent_vector(poly, j))
    end
    exps
end

function extract_exponents(representation::Representation{M}, ring::PolyRing, poly::P) where {M, P<:AbstractAlgebra.Generic.PolyElem}
    exps = Vector{M}(undef, 0)
    @inbounds while !iszero(poly)
        push!(exps, make_ev(M, [AbstractAlgebra.degree(poly)]))
        poly = AbstractAlgebra.tail(poly)
    end
    exps
end

function extract_exponents(representation::Representation{M}, ring::PolyRing, orig_polys::Vector{T}) where {M, T}
    npolys = length(orig_polys)
    exps = Vector{Vector{M}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = extract_exponents(representation, ring, poly)
    end
    exps
end

function exponents_wrt_vars(t, var2idx)
    exp = zeros(Int, length(var2idx))
    @inbounds for (v, p) in Groebner.MultivariatePolynomials.powers(t)
        exp[var2idx[v]] = p
    end
    exp
end

multivariate_length(p::MultivariatePolynomials.AbstractMonomialLike) = 1
multivariate_length(p::MultivariatePolynomials.AbstractTermLike) = 1
multivariate_length(p::AbstractPolynomialLike) = length(p)

function extract_exponents(
            representation::Representation{M},
            ring::PolyRing,
            orig_polys::Vector{T}) where {M, T<:AbstractPolynomialLike{U}} where {U}

    npolys = length(orig_polys)
    exps = Vector{Vector{M}}(undef, npolys)
    vars = MultivariatePolynomials.variables(orig_polys)
    @assert issorted(vars, rev=true)

    var2idx = Dict(vars[i] => i for i in 1:length(vars))
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{M}(undef, multivariate_length(poly))
        @inbounds for (j, t) in enumerate(MultivariatePolynomials.monomials(poly))
            et = exponents_wrt_vars(t, var2idx)
            exps[i][j] = make_ev(M, et)
        end
    end
    exps
end

function extract_exponents(
        representation::Representation{M},
        ring::PolyRing,
        orig_polys::Vector{T},
        ::DegLex) where {M, T}

    npolys = length(orig_polys)
    exps   = Vector{Vector{M}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{M}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = make_ev(M, poly.exps[end-1:-1:1, j])
        end
    end
    exps
end

function extract_exponents(
        representation::Representation{M},
        ring::PolyRing,
        orig_polys::Vector{T},
        ::Lex) where {M, T}

    npolys = length(orig_polys)
    exps   = Vector{Vector{M}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{M}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = make_ev(M, poly.exps[end:-1:1, j])
        end
    end
    exps
end

function extract_exponents(
        representation::Representation{M},
        ring::PolyRing,
        orig_polys::Vector{T},
        ::DegRevLex) where {M, T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{M}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{M}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = make_ev(M, poly.exps[1:end-1, j])
        end
    end
    exps
end

function remove_zeros_from_input!(ring::PolyRing, 
        exps::Vector{Vector{M}}, 
        coeffs::Vector{Vector{T}}) where {M, T}
    @assert length(exps) == length(coeffs)
    filter!(!iszero_coeffvector, coeffs)
    filter!(!iszero_monomvector, exps)
    @assert length(exps) == length(coeffs)
    iszerobasis = isempty(exps)
    if iszerobasis
        push!(exps, Vector{M}())
        push!(coeffs, Vector{T}())
    end
    iszerobasis
end

#------------------------------------------------------------------------------
# Check the consistency of the given monomial ordering
# with respect to the monomial implementation and the length of exponent vector

# Should this be moved to src/monoms ?
@noinline function _throw_monomial_ordering_inconsistent(
        e, o, lo=2, hi=length(e)
    )
    throw(DomainError(o, 
        """The given monomial ordering is inconsistent with the input.
        Exponent: $e
        Indices : $lo to $hi
        Ordering: $o
        Probable cause is that the number of variables does not agree."""
    ))
end

check_ordering(e::M, o::Lex) where {M<:Monom} = true
check_ordering(e::M, o::DegLex) where {M<:Monom} = true
check_ordering(e::M, o::DegRevLex) where {M<:Monom} = true
function check_ordering(e::M, o::Union{Lex, DegLex, DegRevLex}, lo, hi) where {M<:Monom}
    if lo <= hi
        true
    else
        _throw_monomial_ordering_inconsistent(e, wo, lo, hi)
        false
    end
end

function check_ordering(e::M, wo::WeightedOrdering{O}
    ) where {M<:Monom, O<:AbstractMonomialOrdering}
    _throw_monomial_ordering_inconsistent(e, wo)
    false
end
function check_ordering(e::PowerVector{T}, wo::WeightedOrdering{O},
    lo::Int, hi::Int
    ) where {T, O<:AbstractMonomialOrdering}
    check_ordering(e, wo.ord, lo, hi)
    if hi - lo + 1 != length(wo.weights)
        _throw_monomial_ordering_inconsistent(e, wo, lo, hi)
        return false
    end
    true
end
function check_ordering(e::PowerVector{T}, wo::WeightedOrdering{O}
    ) where {T, O<:AbstractMonomialOrdering}
    check_ordering(e, wo, 2, length(e))
end

function check_ordering(e::M, bo::BlockOrdering{R1, R2, O1, O2}
    ) where {M<:Monom, R1, R2, O1<:AbstractMonomialOrdering, O2<:AbstractMonomialOrdering}
    _throw_monomial_ordering_inconsistent(e, bo)
    false
end
function check_ordering(
    e::PowerVector{T}, bo::BlockOrdering{R1, R2, O1, O2},
    lo::Int, hi::Int) where {T, R1, R2, O1<:AbstractMonomialOrdering, O2<:AbstractMonomialOrdering}
    r1 = (first(bo.r1) + 1):(last(bo.r1) + 1)
    r2 = (first(bo.r2) + 1):(last(bo.r2) + 1)
    if first(r1) != lo || last(r2) != hi
        _throw_monomial_ordering_inconsistent(e, bo, lo, hi)
        return false
    end
    check_ordering(e, bo.ord1, first(r1), last(r1))
    check_ordering(e, bo.ord2, first(r2), last(r2))
    true
end
function check_ordering(e::PowerVector{T}, bo::BlockOrdering{R1, R2, O1, O2}
    ) where {T, R1, R2, O1<:AbstractMonomialOrdering, O2<:AbstractMonomialOrdering}
    check_ordering(e, bo, 2, length(e))
end

function check_ordering(e::M, mo::MatrixOrdering) where {M<:Monom}
    _throw_monomial_ordering_inconsistent(e, mo)
    false
end
function check_ordering(e::PowerVector{T}, mo::MatrixOrdering) where {T}
    check_ordering(e, mo, 2, length(e))
end
function check_ordering(
        e::PowerVector{T}, mo::MatrixOrdering,
        lo::Int, hi::Int) where {T}
    rows = mo.rows
    n = hi - lo + 1
    for i in 1:length(rows)
        if length(rows[i]) != n
            _throw_monomial_ordering_inconsistent(e, mo, lo, hi)
            return false
        end
    end
    true
end

#=
    Checks that monomial orderings specified by the given `ring` and `target_ord` 
    are consistent with the given input exponents `exps`.
    In case the target ordering differs from the `ring` ordering,  
    sorts the polynomials terms w.r.t. the target ordering.
    Returns a new polynomial ring in the target ordering.

    Assumes exps and coeffs are non-empty and do not contain zero elements.
=#
function assure_ordering!(
        ring, exps, coeffs, target_ord::O2
    ) where {O2<:AbstractMonomialOrdering}
    @assert !isempty(exps) && !isempty(exps[1])
    check_ordering(exps[1][1], target_ord)
    if ring.ord != target_ord
        sort_input_to_change_ordering!(exps, coeffs, target_ord)
    end
    PolyRing(ring.nvars, target_ord, ring.ch, ring.origring)
end

#------------------------------------------------------------------------------

"""
    Converts internal polynomials for export as elements of `origring`.

    This is the most general implementation.
    It happened to work for polynomials from
        . `Nemo`

    There are also more efficient specializations for:
        . `AbstractAlgebra.Generic.MPoly`
        . `MultivariatePolynomials.AbstractPolynomial`
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {P, M<:Monom, I<:Coeff}

    if ring.origring === :AbstractAlgebra
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs, metainfo)
    elseif ring.origring === :MultivariatePolynomials
        convert_to_output(ring, origpolys, gbexps, gbcoeffs, metainfo)
    elseif ring.origring === :hasparent
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs, metainfo)
    else
        # this actually never happens
        error("This is a bug. Please submit a github issue.")
    end
end

#------------------------------------------------------------------------------

checkexact(c, T::Type{BigInt}) = true
checkexact(c, T::Type{Rational{U}}) where {U} = checkexact(numerator(c), U) && checkexact(denominator(c), U)
function checkexact(c, T)
    if typemin(T) <= c <= typemax(T)
        return true
    else
        throw(DomainError(c, "Coefficient $c in the final basis cannot be converted exactly to $T. Using big arithmetic in input should fix this."))
        return false
    end
end

function check_and_convert_coeffs(coeffs_zz, T)
    cfs = Vector{T}(undef, length(coeffs_zz))
    for i in 1:length(coeffs_zz)
        checkexact(coeffs_zz[i], T)
        cfs[i] = coeffs_zz[i]
    end
    cfs
end

function convert_coeffs_to_output(
        CoeffsVector::Vector{Q},
        ::Type{T}) where {Q<:CoeffQQ, T<:Rational}
    check_and_convert_coeffs(CoeffsVector, T)
end

function convert_coeffs_to_output(
        CoeffsVector::Vector{Q},
        ::Type{T}) where {Q<:CoeffQQ, T<:Integer}
    coeffs_zz = scale_denominators(CoeffsVector)
    check_and_convert_coeffs(coeffs_zz, T)
end

"""
    `multivariate` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {M<:Monom, P<:AbstractPolynomialLike{J}, I<:Coeff} where {J}

    # TODO: hardcoded
    (metainfo.targetord != DegLex()) && @warn "Input polynomial type does not support ordering $(metainfo.targetord). \nComputed basis is correct in $(metainfo.targetord), but terms are ordered in $(DegLex()) in output"

    origvars = MultivariatePolynomials.variables(origpolys)
    # xd
    T = typeof(origpolys[1] + origpolys[1])
    exported = Vector{T}(undef, length(gbexps))
    tmp = Vector{Int}(undef, length(origvars))
    for i in 1:length(gbexps)
        cfs::Vector{J} = convert_coeffs_to_output(gbcoeffs[i], J)
        expvectors = [map(Int, make_dense!(tmp, gbexps[i][j])) for j in 1:length(gbexps[i])]
        expvars = map(t -> t[1]*prod(map(^, origvars, t[2])), zip(cfs, expvectors))
        exported[i] = sum(expvars)
    end
    exported
end

#------------------------------------------------------------------------------

"""
    `hasparent` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origring::R,
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {R, M<:Monom, I}

    @assert hasmethod(base_ring, Tuple{typeof(origring)})

    etype = elem_type(base_ring(origring))
    # rather weak but okay for now
    if etype <: Integer
        coeffs_zz = scale_denominators(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, metainfo)
    else
        convert_to_output(origring, gbexps, gbcoeffs, metainfo)
    end
end

function convert_to_output(
            origring::R,
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {
                R<:Union{AbstractAlgebra.Generic.PolyRing,AbstractAlgebra.PolyRing}, 
                M<:Monom, I
            }

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    @inbounds for i in 1:length(gbexps)
        cfs = zeros(ground, Int(totaldeg(gbexps[i][1]) + 1))
        for (idx, j) in enumerate(gbexps[i])
            cfs[totaldeg(j) + 1] = ground(gbcoeffs[i][idx])
        end
        exported[i] = origring(cfs)
    end
    exported
end

function convert_to_output(
            origring::R,
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {R, M<:Monom, I}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{Int}(undef, AbstractAlgebra.nvars(origring))
    @inbounds for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = [Int.(make_dense!(tmp, gbexps[i][j])) for j in 1:length(gbexps[i])]
        exported[i] = origring(cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------

function create_polynomial(
            origring::AbstractAlgebra.Generic.MPolyRing{T}, coeffs::Vector{T}, exps::Matrix{U}) where {T, U}
    ground = base_ring(origring)
    if !isempty(coeffs)
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring, coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

function create_polynomial(
            origring::AbstractAlgebra.Generic.MPolyRing{T}, coeffs::Vector{T}, exps::Vector{Vector{U}}) where {T, U}
    ground = base_ring(origring)
    if !isempty(coeffs)
        origring(coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

"""
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origring::AbstractAlgebra.Generic.MPolyRing{T},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {M<:Monom,T, I}

    ord = AbstractAlgebra.ordering(origring)
    ordT = ordering_sym2typed(ord)
    if metainfo.targetord != ordT
        ordS = ordering_typed2sym(ord, metainfo.targetord)
        origring, _ = AbstractAlgebra.PolynomialRing(base_ring(origring), AbstractAlgebra.symbols(origring), ordering=ordS)
    end

    if elem_type(base_ring(origring)) <: Integer
        coeffs_zz = scale_denominators(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, metainfo.targetord)
    else
        convert_to_output(origring, gbexps, gbcoeffs, metainfo.targetord)
    end
end

"""
    Rational, Integer, and Finite field degrevlex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::DegRevLex) where  {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{_AA_exponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{_AA_exponenttype}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # TODO x1:
            # Write gbexps[i][jt] directly to exps[1:end-1, jt].
            make_dense!(tmp, gbexps[i][jt])
            exps[1:end-1, jt] .= tmp
            exps[end, jt] = sum(tmp)
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field lex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::Lex) where {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{_AA_exponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{_AA_exponenttype}(undef, nv, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][nv - je + 1]
            # end
            make_dense!(tmp, gbexps[i][jt])
            exps[end:-1:1, jt] .= tmp
        end
        # exps   = UInt64.(hcat(map(x -> x[end-1:-1:1], gbexps[i])...))
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field deglex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::DegLex) where {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{_AA_exponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{_AA_exponenttype}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][nv - je + 1]
            # end
            # exps[nv + 1, jt] = gbexps[i][jt][end]
            make_dense!(tmp, gbexps[i][jt])
            exps[end-1:-1:1, jt] .= tmp
            exps[end, jt] = sum(tmp)
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field, *any other ordering*
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ord::O) where {M, T<:Coeff, U, O<:AbstractMonomialOrdering}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{_AA_exponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Vector{Vector{Int}}(undef, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            make_dense!(tmp, gbexps[i][jt])
            exps[jt] = tmp
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end
