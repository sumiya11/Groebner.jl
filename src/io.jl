
#=
    A note on how we represent polynomials

    First, coefficients, exponents and polynomial ring metainfo
    are extracted from input polynomials with `convert_to_internal`.

    Inside the algorithm all exponents are hashed without collisions,
    so that an integer represents a single monomial.
    That drastically decreases memory consumption.

    For each monomial we also assign a divmask -- a compressed representation
    of exponent vector used to speed up divisibility checks.

    Hence, a polynomial is represented with coefficients vector
    together with vector of hashtable indices.

    Finally, the resulting basis in internal representation is dehashed,
    and passed to `convert_to_output`
=#

#------------------------------------------------------------------------------

"""
    Contains info about polynomial ring
"""
mutable struct PolyRing
    # number of variables
    nvars::Int
    # raw length of exponent vector
    # (should be always equal to nvars+1)
    explen::Int
    # ring monomial ordering,
    # possible are :lex and :degrevlex
    ord::Symbol
    # characteristic of coefficient field
    ch::UInt64
    # information about the original ring of input. Options are:
    #    :abstract for AbstractAlgebra,
    #    :multivariate for MultivariatePolynomials, e.g, DynamicPolynomials,
    #    :hasparent for polynomials constructed with parent ring, e.g., Nemo
    origring::Symbol
end

#------------------------------------------------------------------------------

"""
    Converts input polynomials to internal representation used by the algorithm.
    Extracts base ring information, exponents, and coefficients.

    This is the most general implementation.
    It happened to work for polynomials from
        . `Nemo`

    Currently, there are more efficient specializations for:
        . `AbstractAlgebra.MPoly`
        . `MultivariatePolynomials.AbstractPolynomial`
"""
function convert_to_internal(
            orig_polys::Vector{T},
            ordering::Symbol) where {T}
    isempty(orig_polys) && error("Empty input")
    ordering in (:input, :lex, :degrevlex, :deglex) || error("Not supported ordering $ordering")

    if hasmethod(parent, Tuple{typeof(first(orig_polys))})
        convert_to_internal(orig_polys, ordering, Val(:hasparent))
    else
        error("Sorry, we don't work with this type of polynomials yet. Feel free to open a github issue")
    end
end

function extract_ring(R::T) where {T}
    @assert hasmethod(nvars, Tuple{T})
    # @assert hasmethod(ordering, Tuple{typeof(R)})
    @assert hasmethod(characteristic, Tuple{T})

    nv     = nvars(R)
    explen = nv + 1
    ord = hasmethod(ordering, Tuple{T}) ? ordering(R) : :lex
    ch     = characteristic(R)

    @assert nv + 1 == explen
    @assert ord in (:lex, :degrevlex, :deglex)
    @assert 0 <= ch < 2^32

    PolyRing(nv, explen, ord, UInt64(ch), :hasparent)
end

function extract_coeffs(ring::PolyRing, orig_polys::Vector{T}) where {T}
    if ring.ch > 0
        extract_coeffs_ff(ring, orig_polys)
    else
        extract_coeffs_qq(ring, orig_polys)
    end
end

function extract_coeffs_ff(ring::PolyRing, poly::Poly)
    reverse(map(CoeffFF ∘ data, filter(!iszero, collect(coefficients(poly)))))
end

function extract_coeffs_qq(ring::PolyRing, poly::Poly)
    reverse(map(Rational, filter(!iszero, collect(coefficients(poly)))))
end

function extract_coeffs_ff(ring::PolyRing, poly::MPoly)
    map(CoeffFF ∘ data, coefficients(poly))
end

function extract_coeffs_qq(ring::PolyRing, poly::MPoly)
    map(Rational, coefficients(poly))
end

function extract_coeffs_ff(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{CoeffFF}}(undef, npolys)
    for i in 1:npolys
        coeffs[i] = extract_coeffs_ff(ring, orig_polys[i])
    end
    coeffs
end

function extract_coeffs_qq(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{CoeffQQ}}(undef, npolys)
    for i in 1:npolys
        coeffs[i] = extract_coeffs_qq(ring, orig_polys[i])
    end
    coeffs
end

function extract_exponents(ring, poly::MPoly)
    exps = Vector{ExponentVector}(undef, length(poly))
    @inbounds for j in 1:length(poly)
        exps[j] = ExponentVector(undef, ring.explen)
        exps[j][1:ring.nvars] .= exponent_vector(poly, j)
        exps[j][end] = sum(exps[i][j][k] for k in 1:ring.nvars)
    end
    exps
end

function extract_exponents(ring, poly::Poly)
    # Why define length of univeriate polynomial to be the dense length??
    exps = Vector{ExponentVector}(undef, 0)
    @inbounds while !iszero(poly)
        push!(exps, ExponentVector(undef, ring.explen))
        exps[end][1] = degree(poly)
        exps[end][end] = degree(poly)
        poly = AbstractAlgebra.tail(poly)
    end
    exps
end

function extract_exponents(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    exps = Vector{Vector{ExponentVector}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = extract_exponents(ring, poly)
    end
    exps
end

function extract_polys(ring::PolyRing, orig_polys::Vector{T}) where {T}
    cfs = extract_coeffs(ring, orig_polys)
    exps = extract_exponents(ring, orig_polys)
    exps, cfs
end

function convert_to_internal(
            orig_polys::Vector{T},
            ordering::Symbol,
            ::Val{:hasparent}) where {T}

    R = parent(first(orig_polys))
    ring = extract_ring(R)
    exps, cfs = extract_polys(ring, orig_polys)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

function extract_ring(
        orig_polys::Vector{<:AbstractPolynomialLike{T}}) where {T}

    f = first(orig_polys)

    nv = Groebner.MultivariatePolynomials.nvariables(orig_polys)
    explen = nv + 1
    ord    = :deglex
    ch     = T <: Union{Integer, Rational} ? 0 : characteristic(parent(first(coefficients(f))))

    @assert nv + 1 == explen
    @assert ord in (:lex, :degrevlex, :deglex)
    @assert 0 <= ch < 2^31

    PolyRing(nv, explen, ord, UInt64(ch), :multivariate)
end

function extract_coeffs_qq(
            ring::PolyRing,
            orig_polys::Vector{T}) where {T<:AbstractPolynomialLike{U}} where {U}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{CoeffQQ}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = map(Rational, MultivariatePolynomials.coefficients(poly))
    end
    coeffs
end

function exponents_wrt_vars(t, var2idx)
    exp = zeros(UInt16, length(var2idx))
    @inbounds for (v, p) in Groebner.MultivariatePolynomials.powers(t)
        exp[var2idx[v]] = p
    end
    exp
end

multivariate_length(p::MultivariatePolynomials.AbstractMonomialLike) = 1
multivariate_length(p::MultivariatePolynomials.AbstractTermLike) = 1
multivariate_length(p::AbstractPolynomialLike) = length(p)

function extract_exponents(
            ring::PolyRing,
            orig_polys::Vector{T}) where {T<:AbstractPolynomialLike{U}} where {U}

    npolys = length(orig_polys)
    exps = Vector{Vector{ExponentVector}}(undef, npolys)
    vars = MultivariatePolynomials.variables(orig_polys)
    @assert issorted(vars, rev=true)

    var2idx = Dict(vars[i] => i for i in 1:length(vars))
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{ExponentVector}(undef, multivariate_length(poly))
        @inbounds for (j, t) in enumerate(MultivariatePolynomials.monomials(poly))
            exps[i][j] = ExponentVector(undef, ring.explen)
            exps[i][j][end] = zero(exps[i][j][end])
            et = exponents_wrt_vars(t, var2idx)
            @inbounds for ei in 1:ring.nvars
                exps[i][j][ei] = et[ei]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
    end
    exps
end

"""
    `MultivariatePolynomials.AbstractPolynomial` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{T},
        ordering::Symbol) where {T<:AbstractPolynomialLike{U}} where {U}

    isempty(orig_polys) && error("Empty input")
    ordering in (:input, :lex, :degrevlex, :deglex) || error("Not supported ordering $ordering")

    ring = extract_ring(orig_polys)
    exps, cfs = extract_polys(ring, orig_polys)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

function extract_exponents(
        ring::PolyRing,
        orig_polys::Vector{T},
        ::Ord) where {Ord<:Union{Val{:lex}, Val{:deglex}}} where {T}

    npolys = length(orig_polys)
    exps   = Vector{Vector{ExponentVector}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{ExponentVector}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = ExponentVector(undef, ring.explen)
            exps[i][j][end] = zero(exps[i][j][end])
            @inbounds for ei in 1:ring.nvars
                exps[i][j][ei] = poly.exps[ring.nvars - ei + 1, j]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
    end
    exps
end

function extract_exponents(
        ring::PolyRing,
        orig_polys::Vector{T},
        ::Val{:degrevlex}) where {T}

    npolys = length(orig_polys)
    exps   = Vector{Vector{ExponentVector}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{ExponentVector}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = poly.exps[:, j]
        end
    end
    exps
end

function extract_polys(
        ring::PolyRing,
        orig_polys::Vector{T},
        ord::Symbol) where {T}
    cfs = extract_coeffs(ring, orig_polys)
    exps = extract_exponents(ring, orig_polys, Val(ord))
    exps, cfs
end

"""
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{MPoly{T}},
        ordering::Symbol) where {T}
    isempty(orig_polys) && error("Empty input")

    R = parent(first(orig_polys))
    ring = extract_ring(R)
    exps, cfs = extract_polys(ring, orig_polys, ring.ord)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

"""
    Converts internal polynomials for export as elements of `origring`.

    This is the most general implementation.
    It happened to work for polynomials from
        . `Nemo`

    There are also more efficient specializations for:
        . `AbstractAlgebra.MPoly`
        . `MultivariatePolynomials.AbstractPolynomial`
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {P, I<:Coeff}

    if ring.origring == :abstract
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs, metainfo)
    elseif ring.origring == :multivariate
        convert_to_output(ring, origpolys, gbexps, gbcoeffs, metainfo)
    elseif ring.origring == :hasparent
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs, metainfo)
    else
        # this actually never happens
        error("Sorry, unknown polynomial ring.")
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
        polycoeffs::Vector{CoeffQQ},
        ::Type{T}) where {T<:Rational}
    check_and_convert_coeffs(polycoeffs, T)
end

function convert_coeffs_to_output(
        polycoeffs::Vector{CoeffQQ},
        ::Type{T}) where {T<:Integer}
    coeffs_zz = scale_denominators(polycoeffs)
    check_and_convert_coeffs(coeffs_zz, T)
end

"""
    `multivariate` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {P<:AbstractPolynomialLike{J}, I<:Coeff} where {J}

    # TODO: hardcoded
    (metainfo.targetord != :deglex) && @warn "Input polynomial type does not support ordering $(metainfo.targetord). Computed basis is in $(metainfo.targetord), but terms are ordered in $(:deglex) in output"

    origvars = MultivariatePolynomials.variables(origpolys)
    # xd
    T = typeof(origpolys[1] + origpolys[1])
    exported = Vector{T}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs::Vector{J} = convert_coeffs_to_output(gbcoeffs[i], J)
        expvectors = [map(Int, gbexps[i][j][1:end-1]) for j in 1:length(gbexps[i])]
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
            origring::M,
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {M, T, I}

    @assert hasmethod(base_ring, Tuple{typeof(origring)})

    # TODO: hardcoded
    if hasmethod(ordering, Tuple{M})
        (metainfo.targetord != ordering(origring)) && @warn "Unknown polynomial type. Computed basis is in $(metainfo.targetord), but terms are ordered in $(ordering(origring)) in output"
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
            origring::M,
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {M<:AbstractAlgebra.Generic.PolyRing, I}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    @inbounds for i in 1:length(gbexps)
        cfs = Vector{}
        cfs    = zeros(ground, gbexps[i][1][1] + 1)
        for (idx, j) in enumerate(gbexps[i])
            cfs[j[1] + 1] = ground(gbcoeffs[i][idx])
        end
        exported[i] = origring(cfs)
    end
    exported
end

function convert_to_output(
            origring::M,
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {M<:AbstractAlgebra.Generic.MPolyRing, I}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    @inbounds for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = [Int.(gbexps[i][j][1:end-1]) for j in 1:length(gbexps[i])]
        exported[i] = origring(cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------

function create_polynomial(
            origring::MPolyRing{T}, coeffs::Vector{T}, exps) where {T}
    ground   = base_ring(origring)
    if !isempty(coeffs)
        MPoly{elem_type(ground)}(origring, coeffs, exps)
    else
        MPoly{elem_type(ground)}(origring)
    end
end

"""
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origring::MPolyRing{T},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{I}},
            metainfo::GroebnerMetainfo) where {T, I}

    ord = ordering(origring)
    if metainfo.targetord != ord
        origring, _ = PolynomialRing(base_ring(origring), symbols(origring), ordering=metainfo.targetord)
    end

    if elem_type(base_ring(origring)) <: Integer
        coeffs_zz = scale_denominators(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, Val(ordering(origring)))
    else
        convert_to_output(origring, gbexps, gbcoeffs, Val(ordering(origring)))
    end
end

"""
    Rational, Integer, and Finite field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{U},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:degrevlex}) where  {T<:Coeff, U}

    nv = nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{UInt64}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            for je in 1:nv
                exps[je, jt] = gbexps[i][jt][je]
            end
            exps[nv + 1, jt] = gbexps[i][jt][end]
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field :lex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{U},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:lex}) where {T<:Coeff, U}

    nv = nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{UInt64}(undef, nv, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            for je in 1:nv
                exps[je, jt] = gbexps[i][jt][nv - je + 1]
            end
        end
        # exps   = UInt64.(hcat(map(x -> x[end-1:-1:1], gbexps[i])...))
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

"""
    Rational, Integer, and Finite field :deglex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{U},
            gbexps::Vector{Vector{ExponentVector}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:deglex}) where {T<:Coeff, U}

    nv = nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = Matrix{UInt64}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            for je in 1:nv
                exps[je, jt] = gbexps[i][jt][nv - je + 1]
            end
            exps[nv + 1, jt] = gbexps[i][jt][end]
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------
