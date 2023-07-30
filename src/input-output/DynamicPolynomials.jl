
function peek_at_polynomials(polynomials::Vector{<:AbstractPolynomialLike{T}}) where {T}
    @assert !isempty(polynomials) "Input must not be empty"
    nv = Groebner.MultivariatePolynomials.nvariables(polynomials)
    :dynamicpolynomials, length(polynomials), UInt(0), nv, :deglex
end

function extract_ring(orig_polys::Vector{<:AbstractPolynomialLike{T}}) where {T}
    nv = Groebner.MultivariatePolynomials.nvariables(orig_polys)
    ord = DegLex()
    PolyRing{typeof(ord)}(nv, ord, UInt(0))
end

function extract_coeffs_qq(
    representation,
    ring::PolyRing,
    poly::T
) where {T <: AbstractPolynomialLike{U}} where {U}
    iszero(poly) && (return zero_coeffs_ff(representation.coefftype, ring))
    map(Rational, MultivariatePolynomials.coefficients(poly))
end

function extract_coeffs_qq(
    ring::PolyRing,
    orig_polys::Vector{T}
) where {T <: AbstractPolynomialLike{U}} where {U}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = extract_coeffs_qq(ring, poly)
    end
    coeffs
end

function extract_coeffs(representation, ring, orig_polys)
    extract_coeffs_qq(ring, orig_polys)
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

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T}
) where {T <: AbstractPolynomialLike{U}} where {U}
    M = representation.monomtype
    reversed_order = false
    npolys = length(orig_polys)
    exps = Vector{Vector{M}}(undef, npolys)
    vars = MultivariatePolynomials.variables(orig_polys)
    @assert issorted(vars, rev=true)

    var2idx = Dict(vars[i] => i for i in 1:length(vars))
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{M}(undef, multivariate_length(poly))
        monoms = MultivariatePolynomials.monomials(poly)
        if length(monoms) > 1
            if monoms[1] < monoms[2]
                reversed_order = true
            end
        end
        @inbounds for (j, t) in enumerate(monoms)
            et = exponents_wrt_vars(t, var2idx)
            exps[i][j] = construct_monom(M, et)
        end
    end
    reversed_order, var2idx, exps
end

function convert_to_output(
    ring::PolyRing,
    origpolys::Vector{P},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{I}},
    params::AlgorithmParameters
) where {M <: Monom, P <: AbstractPolynomialLike{J}, I <: Coeff} where {J}
    if isempty(gbexps)
        push!(gbexps, Vector{M}())
        push!(gbcoeffs, Vector{I}())
    end
    if params.target_ord != DegLex()
        @log level = 1 """
        Input polynomial type does not support ordering $(params.target_ord).
        Computed basis is correct in ordering $(params.target_ord), but terms are ordered in DegLex in output"""
    end

    origvars = MultivariatePolynomials.variables(origpolys)
    # xd
    T = typeof(origpolys[1] + origpolys[1])
    exported = Vector{T}(undef, length(gbexps))
    tmp = Vector{Int}(undef, length(origvars))
    for i in 1:length(gbexps)
        cfs::Vector{J} = convert_coeffs_to_output(gbcoeffs[i], J)
        expvectors = [
            map(Int, monom_to_dense_vector!(tmp, gbexps[i][j])) for j in 1:length(gbexps[i])
        ]
        expvars = map(t -> t[1] * prod(map(^, origvars, t[2])), zip(cfs, expvectors))
        exported[i] = sum(expvars)
    end
    exported
end
