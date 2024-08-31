# This file is a part of Groebner.jl. License is GNU GPL v2.

module GroebnerDynamicPolynomialsExt

using DynamicPolynomials

if isdefined(Base, :get_extension)
    using Groebner
else
    using ..Groebner
end

###
# Conversion from DynamicPolynomials.jl to intermediate representation and back.

function Groebner.io_convert_polynomials_to_ir(polynomials::Vector{<:AbstractPolynomialLike{T}}, options::Groebner.KeywordArguments) where {T}
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring = io_extract_ring(polynomials)
    coeffs = io_extract_coeffs(ring, polynomials)
    reversed_order, var_to_index, monoms = io_extract_monoms(ring, polynomials)
    if reversed_order
        for i in 1:length(monoms)
            reverse!(monoms[i])
            reverse!(coeffs[i])
        end
    end
    ring = Groebner.PolyRing(ring.nvars, Groebner.ordering_transform(ring.ord, var_to_index), ring.ch)
    options.ordering = Groebner.ordering_transform(options.ordering, var_to_index)
    ring, monoms, coeffs, options
end

function Groebner.io_convert_ir_to_polynomials(
    ring::Groebner.PolyRing,
    polynomials::Vector{<:AbstractPolynomialLike{T}},
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    options
) where {M <: Groebner.Monom, C <: Groebner.Coeff, T}
    _io_convert_to_output(ring, polynomials, monoms, coeffs, options)
end

function dp_ordering_sym2typed(ord::Symbol)
    if !(ord in (:lex, :deglex, :degrevlex))
        throw(DomainError("Not a supported ordering."))
    end
    if ord === :lex
        Lex()
    elseif ord === :deglex
        DegLex()
    elseif ord === :degrevlex
        DegRevLex()
    end
end

function dp_ord_to_symbol(ord)
    if ord === MultivariatePolynomials.LexOrder
        :lex
    elseif ord === MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}
        :deglex
    elseif ord === MultivariatePolynomials.Graded{
        MultivariatePolynomials.Reverse{MultivariatePolynomials.InverseLexOrder}
    }
        :degrevlex
    else
        throw(DomainError(ord,
            """This ordering from DynamicPolynomials.jl is not supported by Groebner.jl. 
            Consider opening a Github issue if you need it.""")
        )
    end
end

function io_extract_ring(orig_polys::Vector{<:AbstractPolynomialLike{T}}) where {T}
    nv = MultivariatePolynomials.nvariables(orig_polys)
    ord = dp_ord_to_symbol(MultivariatePolynomials.ordering(orig_polys[1]))
    ord_typed = dp_ordering_sym2typed(ord)
    Groebner.PolyRing{typeof(ord_typed), UInt}(nv, ord_typed, UInt(0))
end

function io_extract_coeffs(ring, orig_polys)
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = map(Rational, MultivariatePolynomials.coefficients(poly))
    end
    coeffs
end

function exponents_wrt_vars(t, var2idx)
    exp = zeros(Int, length(var2idx))
    @inbounds for (v, p) in MultivariatePolynomials.powers(t)
        exp[var2idx[v]] = p
    end
    exp
end

multivariate_length(p::MultivariatePolynomials.AbstractMonomialLike) = 1
multivariate_length(p::MultivariatePolynomials.AbstractTermLike) = 1
multivariate_length(p::AbstractPolynomialLike) = MultivariatePolynomials.nterms(p)

function io_extract_monoms(
    ring::Groebner.PolyRing,
    orig_polys::Vector{T}
) where {T <: AbstractPolynomialLike{U}} where {U}
    reversed_order = false
    npolys = length(orig_polys)
    exps = Vector{Vector{Vector{Int}}}(undef, npolys)
    vars = MultivariatePolynomials.variables(orig_polys)
    @assert issorted(vars, rev=true)

    var2idx = Dict(vars[i] => i for i in 1:length(vars))
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{Int}}(undef, multivariate_length(poly))
        monoms = MultivariatePolynomials.monomials(poly)
        if length(monoms) > 1
            if monoms[1] < monoms[2]
                reversed_order = true
            end
        end
        @inbounds for (j, t) in enumerate(monoms)
            exps[i][j] = exponents_wrt_vars(t, var2idx)
        end
    end
    reversed_order, var2idx, exps
end

function _io_convert_to_output(
    ring::Groebner.PolyRing,
    origpolys::Vector{P},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{I}},
    params
) where {M <: Groebner.Monom, P <: AbstractPolynomialLike{J}, I <: Groebner.Coeff} where {J}
    origvars = MultivariatePolynomials.variables(origpolys)
    exported = Vector{Any}(undef, length(gbexps))
    tmp = Vector{Int}(undef, length(origvars))
    for i in 1:length(gbexps)
        cfs = Vector{Rational{BigInt}}(undef, length(gbcoeffs[i]))
        for j in 1:length(gbcoeffs[i])
            cfs[j] = gbcoeffs[i][j]
        end
        expvectors =
            [gbexps[i][j] for j in 1:length(gbexps[i])]
        expvars = map(t -> t[1] * prod(map(^, origvars, t[2])), zip(cfs, expvectors))
        exported[i] = sum(expvars; init=zero(origpolys[1]))
    end
    map(identity, exported)
end


end
