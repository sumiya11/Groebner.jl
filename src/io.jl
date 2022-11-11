
#=
    Our conventions:
    - Trying to compute a Groebner basis of an empty set is an error.
    - The Groebner basis of [0] is [0].
    - The Groebner basis of [0,..., 0] is [0].
    - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]
=#

#=
    A note on how we represent polynomials

    First, coefficients, exponents and polynomial ring
    are extracted from input polynomials with `convert_to_internal`.

    Inside the algorithm all exponents are hashed without collisions,
    so that an integer represents a single monomial.

    Hence, a polynomial is represented with coefficients vector
    together with vector of hashtable indices.

    After the basis is computed, hash table indices are converted back to 
    exponent vectors, and the `convert_to_output` is called.
=#

# The type of exponent vector entries used internally in
# AbstractAlgebra.jl
const AAexponenttype = UInt64

#------------------------------------------------------------------------------

const _supported_orderings = (:lex, :deglex, :degrevlex)

"""
    Contains info about polynomial ring
"""
mutable struct PolyRing{Ch<:CoeffFF}
    # number of variables
    nvars::Int
    # ring monomial ordering,
    # possible are :lex and :degrevlex
    ord::Symbol
    # characteristic of coefficient field
    ch::Ch
    # information about the original ring of input. Options are:
    #    :AbstractAlgebra for AbstractAlgebra,
    #    :MultivariatePolynomials for MultivariatePolynomials, e.g, DynamicPolynomials,
    #    :hasparent for polynomials constructed with parent ring, e.g., Nemo
    #    :undefined
    origring::Symbol
end

#------------------------------------------------------------------------------

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

"""
    Converts input polynomials to internal representation used by the algorithm.
    Extracts base ring information, exponents, and coefficients.

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
        ordering::Symbol) where {T}
    isempty(orig_polys) && throw(DomainError(orig_polys, "Empty input."))
    ordering in _supported_orderings || ordering === :input || throw(DomainError(ordering, "Not supported ordering."))

    if hasmethod(AbstractAlgebra.parent, Tuple{typeof(first(orig_polys))})
        convert_to_internal(representation, orig_polys, ordering, Val(:hasparent))
    else
        convert_to_internal(representation, orig_polys, ordering, Val(:undefined))
    end
end

function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::Symbol,
        ::Val{:undefined}) where {T}
    error("Sorry, we don't work with this type of polynomials ($T) yet. Feel free to open an issue")
end

"""
    `hasparent` convention specialization
"""
function convert_to_internal(
        representation,
        orig_polys::Vector{T},
        ordering::Symbol,
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
        ordering::Symbol,
        ::Val{:undefined}) where {T<:AbstractPolynomialLike{U}} where {U}
    ring = extract_ring(orig_polys)
    check_domain(representation, ring)
    exps, cfs = extract_polys(representation, ring, orig_polys, ring.ord)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

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
    @assert ord in _supported_orderings
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
    ord = hasmethod(AbstractAlgebra.ordering, Tuple{T}) ? AbstractAlgebra.ordering(R) : :deglex
    ch     = AbstractAlgebra.characteristic(R)

    check_ground_domain(nv, ord, ch)

    Ch = determinechartype(ch)

    PolyRing{Ch}(nv, ord, Ch(ch), :AbstractAlgebra)
end

function extract_ring(orig_polys::Vector{<:AbstractPolynomialLike{T}}) where {T}
    f = first(orig_polys)

    nv = Groebner.MultivariatePolynomials.nvariables(orig_polys)
    ord    = :deglex
    ch     = 0

    check_ground_domain(nv, ord, ch)

    Ch = determinechartype(ch)

    PolyRing{Ch}(nv, ord, Ch(ch), :MultivariatePolynomials)
end

#------------------------------------------------------------------------------

function extract_polys(
        representation,
        ring::PolyRing,
        orig_polys::Vector{T},
        ord::Symbol) where {T}
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
function extract_coeffs_ff(ring::PolyRing{Ch}, poly::AbstractAlgebra.Generic.Poly) where {Ch}
    iszero(poly) && (return zero_coeffvector_ff(ring))
    reverse(map(Ch ∘ AbstractAlgebra.data, filter(!iszero, collect(AbstractAlgebra.coefficients(poly)))))
end

# specialization for univariate polynomials
function extract_coeffs_qq(ring::PolyRing, poly::AbstractAlgebra.Generic.Poly)
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
        ::Val{:deglex}) where {M, T}

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
        ::Val{:lex}) where {M, T}

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
        ::Val{:degrevlex}) where {M, T}
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

function assure_ordering!(ring, exps, coeffs, target_ord)
    if ring.ord != target_ord
        sort_input_to_change_ordering!(exps, coeffs, target_ord)
    end
    ring.ord = target_ord
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

    if ring.origring == :AbstractAlgebra
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs, metainfo)
    elseif ring.origring == :MultivariatePolynomials
        convert_to_output(ring, origpolys, gbexps, gbcoeffs, metainfo)
    elseif ring.origring == :hasparent
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
    (metainfo.targetord != :deglex) && @warn "Input polynomial type does not support ordering $(metainfo.targetord). \nComputed basis is correct in $(metainfo.targetord), but terms are ordered in $(:deglex) in output"

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

    # TODO: hardcoded
    if hasmethod(AbstractAlgebra.ordering, Tuple{R})
        (metainfo.targetord != AbstractAlgebra.ordering(origring)) && @warn "Unknown polynomial type. Computed basis is in $(metainfo.targetord), but terms are ordered in $(ordering(origring)) in output"
    end

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
            metainfo::GroebnerMetainfo) where {R<:AbstractAlgebra.Generic.PolyRing, M<:Monom, I}

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
            origring::AbstractAlgebra.Generic.MPolyRing{T}, coeffs::Vector{T}, exps) where {T}
    ground = base_ring(origring)
    if !isempty(coeffs)
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring, coeffs, exps)
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
    if metainfo.targetord != ord
        origring, _ = AbstractAlgebra.PolynomialRing(base_ring(origring), AbstractAlgebra.symbols(origring), ordering=metainfo.targetord)
    end

    if elem_type(base_ring(origring)) <: Integer
        coeffs_zz = scale_denominators(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, Val(AbstractAlgebra.ordering(origring)))
    else
        convert_to_output(origring, gbexps, gbcoeffs, Val(AbstractAlgebra.ordering(origring)))
    end
end

"""
    Rational, Integer, and Finite field :degrevlex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:degrevlex}) where  {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{AAexponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{AAexponenttype}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][je]
            # end
            # exps[nv + 1, jt] = gbexps[i][jt][end]
            make_dense!(tmp, gbexps[i][jt])
            exps[1:end-1, jt] .= tmp
            exps[end, jt] = sum(tmp)
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field :lex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:lex}) where {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{AAexponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{AAexponenttype}(undef, nv, length(gbcoeffs[i]))
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
    Rational, Integer, and Finite field :deglex
    `AbstractAlgebra.Generic.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::AbstractAlgebra.Generic.MPolyRing{U},
            gbexps::Vector{Vector{M}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:deglex}) where {M, T<:Coeff, U}

    nv = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp = Vector{AAexponenttype}(undef, nv)
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{AAexponenttype}(undef, nv + 1, length(gbcoeffs[i]))
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

#------------------------------------------------------------------------------
