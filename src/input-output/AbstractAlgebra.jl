# Conversion from AbstractAlgebra.jl into internal representation.

# This all is just not fantastic. 
# We are practically dancing for rain around AbstractAlgebra.jl internals here.

const _AA_supported_orderings_symbols = (:lex, :deglex, :degrevlex)
const _AA_exponent_type = UInt64

function peek_at_polynomials(polynomials::Vector{T}) where {T}
    @assert !isempty(polynomials) "Input must not be empty"
    R = parent(first(polynomials))
    nvars = AbstractAlgebra.nvars(R)
    ord = if hasmethod(AbstractAlgebra.ordering, Tuple{typeof(R)})
        # if multivariate
        AbstractAlgebra.ordering(R)
    else
        # if univariate, defaults to lex
        :lex
    end
    char = Int(AbstractAlgebra.characteristic(R))
    :abstractalgebra, length(polynomials), char, nvars, ord
end

# Determines the monomial ordering of the output,
# given the original ordering `origord` and the targer ordering `targetord`
ordering_typed2sym(origord, targetord::Lex) = :lex
ordering_typed2sym(origord, targetord::DegLex) = :deglex
ordering_typed2sym(origord, targetord::DegRevLex) = :degrevlex
ordering_typed2sym(origord) = origord
ordering_typed2sym(origord, targetord::AbstractMonomialOrdering) = origord

function ordering_sym2typed(ord::Symbol)
    ord in _AA_supported_orderings_symbols ||
        throw(DomainError(ord, "Not a supported ordering."))
    if ord === :lex
        Lex()
    elseif ord === :deglex
        DegLex()
    elseif ord === :degrevlex
        DegRevLex()
    end
end

function extract_ring(polynomials)
    R = parent(first(polynomials))
    T = typeof(R)
    @assert hasmethod(AbstractAlgebra.nvars, Tuple{T})
    @assert hasmethod(AbstractAlgebra.characteristic, Tuple{T})
    nv = AbstractAlgebra.nvars(R)
    # lex is the default ordering on univariate polynomials
    ord =
        hasmethod(AbstractAlgebra.ordering, Tuple{T}) ? AbstractAlgebra.ordering(R) :
        :lex
    # type unstable:
    ordT = ordering_sym2typed(ord)
    ch   = AbstractAlgebra.characteristic(R)
    PolyRing{typeof(ordT)}(nv, ordT, Int(BigInt(ch)))
end

function extract_coeffs(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T}
) where {T}
    npolys = length(orig_polys)
    coeffs = Vector{Vector{representation.coefftype}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        coeffs[i] = extract_coeffs(representation, ring, poly)
    end
    coeffs
end

function extract_coeffs(representation::PolynomialRepresentation, ring::PolyRing, orig_poly)
    if ring.ch > 0
        extract_coeffs_ff(representation, ring, orig_poly)
    else
        extract_coeffs_qq(representation, ring, orig_poly)
    end
end

# specialization for univariate polynomials
function extract_coeffs_ff(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    poly::Union{AbstractAlgebra.Generic.Poly, AbstractAlgebra.PolyElem}
)
    iszero(poly) && (return zero_coeffs_ff(representation.coefftype, ring))
    reverse(
        map(
            # NOTE: type instable!
            representation.coefftype ∘ AbstractAlgebra.data,
            filter(!iszero, collect(AbstractAlgebra.coefficients(poly)))
        )
    )
end

# specialization for univariate polynomials
function extract_coeffs_qq(
    representation,
    ring::PolyRing,
    poly::Union{AbstractAlgebra.Generic.Poly, AbstractAlgebra.PolyElem}
)
    iszero(poly) && (return zero_coeffs_qq(representation.coefftype, ring))
    reverse(map(Rational, filter(!iszero, collect(AbstractAlgebra.coefficients(poly)))))
end

# specialization for multivariate polynomials
function extract_coeffs_ff(
    representation::PolynomialRepresentation,
    ring::PolyRing{T},
    poly
) where {T}
    iszero(poly) && (return zero_coeffs_ff(representation.coefftype, ring))
    Ch = representation.coefftype
    # TODO: Get rid of this composed function !
    map(Ch ∘ AbstractAlgebra.data, AbstractAlgebra.coefficients(poly))
end

function extract_coeffs_raw!(graph, representation, polys::Vector{T}, kws) where {T}
    # write directly to graph.basis
    # TODO: 2dim permutation to handle different orderings!
    ring = extract_ring(polys)
    basis = graph.buf_basis
    perm = graph.input_permutation
    Ch = representation.coefftype
    _extract_coeffs_raw!(basis, perm, polys, Ch)
    @log level=-5 "Extracted coefficients from $(length(polys)) polynomials." basis perm
    ring
end

function _extract_coeffs_raw!(basis, perm, polys, ::Type{CoeffsType}) where {CoeffsType}
    @assert basis.nfilled == count(!iszero, polys)
    polys = filter(!iszero, polys)
    @inbounds for i in 1:length(polys)
        poly = polys[perm[i]]
        basis_cfs = basis.coeffs[i]
        @assert length(poly) == length(basis_cfs)
        for j in 1:length(poly)
            basis_cfs[j] =
                convert(CoeffsType, AbstractAlgebra.data(AbstractAlgebra.coeff(poly, j)))  # TODO: PERMUTE
        end
    end
    nothing
end

# specialization for multivariate polynomials
function extract_coeffs_qq(representation, ring::PolyRing, poly)
    iszero(poly) && (return zero_coeffs_qq(representation.coefftype, ring))
    map(Rational, AbstractAlgebra.coefficients(poly))
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    poly
)
    exps = Vector{representation.monomtype}(undef, length(poly))
    @inbounds for j in 1:length(poly)
        exps[j] = construct_monom(
            representation.monomtype,
            AbstractAlgebra.exponent_vector(poly, j)
        )
    end
    exps
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    poly::P
) where {P <: AbstractAlgebra.Generic.PolyElem}
    exps = Vector{representation.monomtype}(undef, 0)
    @inbounds while !iszero(poly)
        push!(
            exps,
            construct_monom(representation.monomtype, [AbstractAlgebra.degree(poly)])
        )
        poly = AbstractAlgebra.tail(poly)
    end
    exps
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T}
) where {T}
    npolys = length(orig_polys)
    exps = Vector{Vector{representation.monomtype}}(undef, npolys)
    @inbounds for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = extract_monoms(representation, ring, poly)
    end
    exps
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T},
    ::DegLex
) where {T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{representation.monomtype}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{representation.monomtype}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] =
                construct_monom(representation.monomtype, poly.exps[(end - 1):-1:1, j])
        end
    end
    exps
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T},
    ::Lex
) where {T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{representation.monomtype}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{representation.monomtype}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = construct_monom(representation.monomtype, poly.exps[end:-1:1, j])
        end
    end
    exps
end

function extract_monoms(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    orig_polys::Vector{T},
    ::DegRevLex
) where {T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{representation.monomtype}}(undef, npolys)
    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{representation.monomtype}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] =
                construct_monom(representation.monomtype, poly.exps[1:(end - 1), j])
        end
    end
    exps
end

function convert_to_output(
    origring::R,
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{I}},
    params::AlgorithmParameters
) where {
    R <: Union{AbstractAlgebra.Generic.PolyRing, AbstractAlgebra.PolyRing},
    M <: Monom,
    I
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
    params::AlgorithmParameters
) where {R, M <: Monom, I}
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{Int}(undef, AbstractAlgebra.nvars(origring))
    @inbounds for i in 1:length(gbexps)
        cfs = map(ground, gbcoeffs[i])
        exps =
            [Int.(monom_to_dense_vector!(tmp, gbexps[i][j])) for j in 1:length(gbexps[i])]
        exported[i] = origring(cfs, exps)
    end
    exported
end

function create_polynomial(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    coeffs::Vector{T},
    exps::Matrix{U}
) where {T, U}
    ground = base_ring(origring)
    if !isempty(coeffs)
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring, coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

function create_polynomial(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    coeffs::Vector{T},
    exps::Vector{Vector{U}}
) where {T, U}
    ground = base_ring(origring)
    if !isempty(coeffs)
        origring(coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

function convert_to_output(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{I}},
    params::AlgorithmParameters
) where {M <: Monom, T, I}
    ord = AbstractAlgebra.ordering(origring)
    ordT = ordering_sym2typed(ord)
    if params.target_ord != ordT
        ordS = ordering_typed2sym(ord, params.target_ord)
        origring, _ = AbstractAlgebra.PolynomialRing(
            base_ring(origring),
            AbstractAlgebra.symbols(origring),
            ordering=ordS
        )
    end

    if elem_type(base_ring(origring)) <: Integer
        coeffs_zz = clear_denominators(gbcoeffs)
        convert_to_output(origring, gbexps, coeffs_zz, params.target_ord)
    else
        convert_to_output(origring, gbexps, gbcoeffs, params.target_ord)
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
    ::DegRevLex
) where {M, T <: Coeff, U}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{_AA_exponent_type}(undef, nv)
    for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{_AA_exponent_type}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][je]
            # end
            # exps[nv + 1, jt] = gbexps[i][jt][end]
            monom_to_dense_vector!(tmp, gbexps[i][jt])
            exps[1:(end - 1), jt] .= tmp
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
    ::Lex
) where {M, T <: Coeff, U}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{_AA_exponent_type}(undef, nv)
    for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{_AA_exponent_type}(undef, nv, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][nv - je + 1]
            # end
            monom_to_dense_vector!(tmp, gbexps[i][jt])
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
    ::DegLex
) where {M, T <: Coeff, U}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{_AA_exponent_type}(undef, nv)
    for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{_AA_exponent_type}(undef, nv + 1, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            # for je in 1:nv
            #     exps[je, jt] = gbexps[i][jt][nv - je + 1]
            # end
            # exps[nv + 1, jt] = gbexps[i][jt][end]
            monom_to_dense_vector!(tmp, gbexps[i][jt])
            exps[(end - 1):-1:1, jt] .= tmp
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
    ord::O
) where {M, T <: Coeff, U, O <: AbstractMonomialOrdering}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{_AA_exponent_type}(undef, nv)
    for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Vector{Vector{Int}}(undef, length(gbcoeffs[i]))
        @inbounds for jt in 1:length(gbcoeffs[i])
            monom_to_dense_vector!(tmp, gbexps[i][jt])
            exps[jt] = tmp
        end
        exported[i] = create_polynomial(origring, cfs, exps)
    end
    exported
end
