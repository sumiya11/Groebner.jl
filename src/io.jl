
#=
    A note on how we represent exponents internally.


=#

#------------------------------------------------------------------------------

"""
    Contains info about polynomial ring
"""
mutable struct PolyRing
    # number of variables
    nvars::Int
    # raw length of exponent vector
    explen::Int
    # ring monomial ordering,
    # possible are :lex and :degrevlex
    ord::Symbol
    # characteristic of coefficient field
    ch::UInt64
    # information about the original ring of input. Options are:
    #    :abstract for AbstractAlgebra,
    #    :multivariate for MultivariatePolynomials,
    #    :hasparent for polynomials constructed with parent ring
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
    @assert hasmethod(nvars, Tuple{typeof(R)})
    @assert hasmethod(ordering, Tuple{typeof(R)})
    @assert hasmethod(characteristic, Tuple{typeof(R)})

    nv     = nvars(R)
    explen = nv + 1
    ord    = ordering(R)
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

function extract_coeffs_ff(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)
    for i in 1:npolys
        coeffs[i] = map(UInt64 ∘ data, coefficients(orig_polys[i]))
    end
    coeffs
end

function extract_coeffs_qq(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)
    for i in 1:npolys
        coeffs[i] = map(Rational, coefficients(orig_polys[i]))
    end
    coeffs
end

function extract_exponents(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    exps = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, ring.explen)
            exps[i][j][1:ring.nvars] .= exponent_vector(poly, j)
            exps[i][j][end] = sum(exps[i][j][k] for k in 1:ring.nvars)
        end
    end
    exps
end

function extract_polys(ring::PolyRing, orig_polys::Vector{T}) where {T}
    cfs = extract_coeffs(ring, orig_polys)
    exps = extract_exponents(ring, orig_polys)
    exps, cfs
end

function assure_ordering!(ring, ordering, exps, cfs)
    if ordering != :input
        if ordering != ring.ord
            ring.ord = ordering
            sort_input_to_change_ordering!(ordering, exps, cfs)
        end
    end
end

function convert_to_internal(
            orig_polys::Vector{T},
            ordering::Symbol,
            ::Val{:hasparent}) where {T}

    R = parent(first(orig_polys))
    ring = extract_ring(R)
    exps, cfs = extract_polys(ring, orig_polys)
    assure_ordering!(ring, ordering, exps, cfs)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

# TODO: check type stability
function extract_ring(
        orig_polys::Vector{<:AbstractPolynomial{T}}) where {T}

    f = first(orig_polys)

    nv = maximum(map(Groebner.MultivariatePolynomials.nvariables, orig_polys))
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
            orig_polys::Vector{T}) where {T<:AbstractPolynomial{U}} where {U}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = map(Rational, MultivariatePolynomials.coefficients(poly))
    end
    coeffs
end

function exponents_wrt_vars(t, var2idx)
    exp = zeros(UInt16, length(var2idx))
    for (v, p) in Groebner.MultivariatePolynomials.powers(t)
        exp[var2idx[v]] = p
    end
    exp
end

function extract_exponents(
            ring::PolyRing,
            orig_polys::Vector{T}) where {T<:AbstractPolynomial{U}} where {U}

    npolys = length(orig_polys)
    exps = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    vars = MultivariatePolynomials.variables(orig_polys)
    @assert issorted(vars, rev=true)

    var2idx = Dict(vars[i] => i for i in 1:length(vars))
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for (j, t) in enumerate(Groebner.MultivariatePolynomials.monomials(poly))
            exps[i][j] = Vector{UInt16}(undef, ring.explen)
            exps[i][j][end] = zero(exps[i][j][end])
            et = exponents_wrt_vars(t, var2idx)
            for ei in 1:ring.nvars
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
        ordering::Symbol) where {T<:AbstractPolynomial{U}} where {U}

    isempty(orig_polys) && error("Empty input")
    ordering in (:input, :lex, :degrevlex, :deglex) || error("Not supported ordering $ordering")

    ring = extract_ring(orig_polys)
    exps, cfs = extract_polys(ring, orig_polys)
    assure_ordering!(ring, ordering, exps, cfs)
    ring, exps, cfs
end

#------------------------------------------------------------------------------

function extract_exponents(
        ring::PolyRing,
        orig_polys::Vector{T},
        ::Ord) where {Ord<:Union{Val{:lex}, Val{:deglex}}} where {T}

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, ring.explen)
            exps[i][j][end] = zero(exps[i][j][end])
            for ei in 1:ring.nvars
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
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
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
    assure_ordering!(ring, ordering, exps, cfs)
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
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {P, I}

    if ring.origring == :abstract
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs)
    elseif ring.origring == :multivariate
        convert_to_output(ring, origpolys, gbexps, gbcoeffs)
    elseif ring.origring == :hasparent
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs)
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
        # TODO
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
        polycoeffs::Vector{Rational{BigInt}},
        ::Type{T}) where {T<:Rational}
    check_and_convert_coeffs(polycoeffs, T)
end

function convert_coeffs_to_output(
        polycoeffs::Vector{Rational{BigInt}},
        ::Type{T}) where {T<:Integer}
    coeffs_zz = scale_denominators!(polycoeffs)
    check_and_convert_coeffs(coeffs_zz, T)
end

"""
    `multivariate` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {P<:AbstractPolynomial{J}, I<:Rational} where {J}

    origvars = MultivariatePolynomials.variables(first(origpolys))
    exported = Vector{P}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs::Vector{J} = convert_coeffs_to_output(gbcoeffs[i], J)
        expvectors = [gbexps[i][j][1:end-1] for j in 1:length(gbexps[i])]
        expvars = map(ev -> prod(origvars .^ Int.(ev)), expvectors)
        exported[i] = P(cfs, expvars)
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
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {M, T, I}

    @assert hasmethod(base_ring, Tuple{typeof(origring)})

    etype = elem_type(base_ring(origring))
    # rather weak but okay for now
    if etype <: Integer
        coeffs_zz = scale_denominators!(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz)
    else
        convert_to_output(origring, gbexps, gbcoeffs)
    end
end

function convert_to_output(
            origring::M,
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where{M, I}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
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
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {T, I}

    ord = ordering(origring)
    if elem_type(base_ring(origring)) <: Integer
        coeffs_zz = scale_denominators!(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, Val(ord))
    else
        convert_to_output(origring, gbexps, gbcoeffs, Val(ord))
    end
end

"""
    Rational, Integer, and Finite field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{U},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:degrevlex}) where  {T<:Union{Rational{BigInt},BigInt,UInt64}, U}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(gbexps[i]...))
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
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:lex}) where {T<:Union{Rational{BigInt},BigInt,UInt64}, U}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> x[end-1:-1:1], gbexps[i])...))
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
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{T}},
            ::Val{:deglex}) where {T<:Union{Rational{BigInt},BigInt,UInt64}, U}

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> [x[end-1:-1:1]..., x[end]], gbexps[i])...))
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------
