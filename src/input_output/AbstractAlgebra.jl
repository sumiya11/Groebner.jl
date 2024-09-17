# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Conversion from AbstractAlgebra.jl/Nemo.jl to intermediate representation and back.

# Conversions in this file must be thread-safe, since functions in the interface
# can be used in parallel. In particular, we rely on the fact that some
# functions in AbstractAlgebra.jl are thread-safe
# https://github.com/Nemocas/AbstractAlgebra.jl/issues/1542

const aa_supported_orderings = (:lex, :deglex, :degrevlex)
const aa_exponent_type = UInt64

aa_is_multivariate_ring(ring) = AbstractAlgebra.elem_type(ring) <: AbstractAlgebra.MPolyRingElem

function io_convert_polynomials_to_ir(polynomials, options::KeywordArguments)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring = io_extract_ring(polynomials)
    coeffs = io_extract_coeffs_ir(ring, polynomials)
    reversed_order, var_to_index, monoms = io_extract_monoms_ir(ring, polynomials)
    @invariant length(coeffs) == length(monoms)
    ring = PolyRing(ring.nvars, ordering_transform(ring.ord, var_to_index), ring.ch)
    options.ordering = ordering_transform(options.ordering, var_to_index)
    ring, monoms, coeffs, options
end

function io_convert_ir_to_polynomials(
    ring::PolyRing,
    polynomials,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    options
) where {M <: Monom, C <: Coeff}
    origring = AbstractAlgebra.parent(first(polynomials))
    _io_convert_ir_to_polynomials(origring, monoms, coeffs, options)::typeof(polynomials)
end

function io_extract_ring(polynomials)
    if !all(
        f -> f isa AbstractAlgebra.MPolyRingElem || f isa AbstractAlgebra.PolyRingElem,
        polynomials
    )
        throw(DomainError("Unknown type."))
    end
    _io_check_input(polynomials)
    R = parent(first(polynomials))
    nv = AbstractAlgebra.nvars(R)
    # lex is the default ordering on univariate polynomials
    ord = if aa_is_multivariate_ring(R)
        AbstractAlgebra.internal_ordering(R)
    else
        :lex
    end
    # type unstable:
    ordT = ordering_sym2typed(ord)
    ch   = AbstractAlgebra.characteristic(R)
    PolyRing{typeof(ordT), UInt}(nv, ordT, UInt(ch))
end

function io_extract_coeffs_ir(ring::PolyRing, polys)
    if ring.ch > 0
        io_extract_coeffs_ir_ff(ring, polys)
    else
        io_extract_coeffs_ir_qq(ring, polys)
    end
end

io_lift_coeff_ff(c) = UInt64(AbstractAlgebra.lift(c))
io_lift_coeff_ff(c::AbstractAlgebra.GFElem) = AbstractAlgebra.data(c)
io_lift_coeff_ff(c::Union{Nemo.FqFieldElem, Nemo.fpFieldElem}) =
    UInt64(AbstractAlgebra.lift(Nemo.ZZ, c))

function io_extract_coeffs_ir_ff(ring::PolyRing{T}, polys) where {T}
    res = Vector{Vector{UInt64}}(undef, length(polys))
    for i in 1:length(polys)
        poly = polys[i]
        if !aa_is_multivariate_ring(parent(polys[1]))
            res[i] = map(io_lift_coeff_ff, AbstractAlgebra.coefficients(poly))
            reverse!(res[i])
            filter!(!iszero, res[i])
        else
            res[i] = map(io_lift_coeff_ff, AbstractAlgebra.coefficients(poly))
        end
    end
    res
end

function io_extract_coeffs_ir_qq(ring::PolyRing, polys)
    res = Vector{Vector{Rational{BigInt}}}(undef, length(polys))
    for i in 1:length(polys)
        poly = polys[i]
        if !aa_is_multivariate_ring(parent(polys[1]))
            res[i] = map(Rational{BigInt}, collect(AbstractAlgebra.coefficients(poly)))
            reverse!(res[i])
            filter!(!iszero, res[i])
        else
            res[i] = map(Rational{BigInt}, collect(AbstractAlgebra.coefficients(poly)))
        end
    end
    res
end

function io_extract_monoms_ir(ring::PolyRing, polys)
    ring_aa = AbstractAlgebra.parent(polys[1])
    v = AbstractAlgebra.gens(ring_aa)
    var_to_index = Dict{elem_type(ring_aa), Int}(v .=> 1:AbstractAlgebra.nvars(ring_aa))
    res = Vector{Vector{Vector{UInt64}}}(undef, length(polys))
    @inbounds for i in 1:length(polys)
        poly = polys[i]
        if !aa_is_multivariate_ring(parent(polys[1]))
            perm = filter(
                j -> !iszero(AbstractAlgebra.coeff(polys[i], j - 1)),
                collect(1:length(polys[i]))
            )
            res[i] = collect(AbstractAlgebra.exponent_vectors(polys[i]))
            res[i] = res[i][perm]
            reverse!(res[i])
        else
            res[i] = collect(AbstractAlgebra.exponent_vectors(polys[i]))
        end
    end
    false, var_to_index, res
end

###
# Converting from AbstractAlgebra to internal representation.

function _io_check_input(polynomials::Vector{T}) where {T}
    R = AbstractAlgebra.parent(first(polynomials))
    K = AbstractAlgebra.base_ring(R)
    if !(K isa AbstractAlgebra.Field)
        throw(DomainError("Coefficient ring must be a field: $K"))
    end
    if (AbstractAlgebra.characteristic(K) > typemax(UInt64))
        throw(DomainError("Field characteristic must be less than 2^64"))
    end
    if !iszero(AbstractAlgebra.characteristic(K))
        if !isone(AbstractAlgebra.degree(K))
            throw(DomainError("Non-prime coefficient fields are not supported"))
        end
    end
    true
end

function ordering_sym2typed(ord::Symbol)
    if !(ord in aa_supported_orderings)
        __throw_input_not_supported(ord, "Not a supported ordering.")
    end
    if ord === :lex
        Lex()
    elseif ord === :deglex
        DegLex()
    elseif ord === :degrevlex
        DegRevLex()
    end
end

###
# Converting from internal representation to AbstractAlgebra.jl

# Specialization for univariate polynomials
function _io_convert_ir_to_polynomials(
    origring::R,
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    params
) where {
    R <: Union{AbstractAlgebra.Generic.PolyRing, AbstractAlgebra.PolyRing},
    M <: Monom,
    C <: Coeff
}
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    @inbounds for i in 1:length(gbexps)
        if isempty(gbexps[i])
            exported[i] = origring()
            continue
        end
        cfs = zeros(ground, Int(sum(gbexps[i][1]) + 1))
        for (idx, j) in enumerate(gbexps[i])
            cfs[sum(j) + 1] = ground(gbcoeffs[i][idx])
        end
        exported[i] = origring(cfs)
    end
    exported
end

# The most generic specialization
# (Nemo.jl opts for this specialization)
function _io_convert_ir_to_polynomials(
    origring::R,
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    params
) where {R, M <: Monom, C <: Coeff}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    @inbounds for i in 1:length(gbexps)
        cfs = map(ground, gbcoeffs[i])
        exps = Vector{Vector{Int}}(undef, length(gbcoeffs[i]))
        for jt in 1:length(gbcoeffs[i])
            exps[jt] = gbexps[i][jt]
        end
        exported[i] = origring(cfs, exps)
    end
    exported
end

function create_aa_polynomial(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    coeffs::Vector{T},
    exps::Matrix{U}
) where {T, U}
    ground = AbstractAlgebra.base_ring(origring)
    if !isempty(coeffs)
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring, coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

function create_aa_polynomial(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    coeffs::Vector{T},
    exps::Vector{Vector{U}}
) where {T, U}
    ground = AbstractAlgebra.base_ring(origring)
    if !isempty(coeffs)
        origring(coeffs, exps)
    else
        AbstractAlgebra.Generic.MPoly{elem_type(ground)}(origring)
    end
end

# Dispatch between monomial orderings
function _io_convert_ir_to_polynomials(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    params
) where {T, M <: Monom, C <: Coeff}
    ord_aa = AbstractAlgebra.internal_ordering(origring)
    _ord_aa = ordering_sym2typed(ord_aa)
    input_ordering_matches_output = true
    target_ord = params.ordering isa InputOrdering ? _ord_aa : params.ordering
    if (target_ord != _ord_aa)
        input_ordering_matches_output = false
    end
    _io_convert_ir_to_polynomials(
        origring,
        gbexps,
        gbcoeffs,
        target_ord,
        Val{input_ordering_matches_output}()
    )
end

# Specialization for degrevlex for matching orderings
function _io_convert_ir_to_polynomials(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    ::DegRevLex,
    input_ordering_matches_output::Val{true}
) where {T, M <: Monom, C <: Coeff}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{aa_exponent_type}(undef, nv)
    @inbounds for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{aa_exponent_type}(undef, nv + 1, length(gbcoeffs[i]))
        for jt in 1:length(gbcoeffs[i])
            for k in 1:length(tmp)
                exps[k, jt] = gbexps[i][jt][k]
            end
            exps[end, jt] = sum(gbexps[i][jt])
        end
        exported[i] = create_aa_polynomial(origring, cfs, exps)
    end
    exported
end

# Specialization for lex for matching orderings
function _io_convert_ir_to_polynomials(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    ::Lex,
    input_ordering_matches_output::Val{true}
) where {T, M <: Monom, C <: Coeff}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{aa_exponent_type}(undef, nv)
    @inbounds for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{aa_exponent_type}(undef, nv, length(gbcoeffs[i]))
        for jt in 1:length(gbcoeffs[i])
            exps[end:-1:1, jt] .= gbexps[i][jt]
        end
        exported[i] = create_aa_polynomial(origring, cfs, exps)
    end
    exported
end

# Specialization for deglex for matching orderings
function _io_convert_ir_to_polynomials(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    ::DegLex,
    input_ordering_matches_output::Val{true}
) where {T, M <: Monom, C <: Coeff}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{aa_exponent_type}(undef, nv)
    @inbounds for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Matrix{aa_exponent_type}(undef, nv + 1, length(gbcoeffs[i]))
        for jt in 1:length(gbcoeffs[i])
            # monom_to_vector!(tmp, gbexps[i][jt])
            exps[(end - 1):-1:1, jt] .= gbexps[i][jt]
            exps[end, jt] = sum(gbexps[i][jt])
        end
        exported[i] = create_aa_polynomial(origring, cfs, exps)
    end
    exported
end

# All other orderings
function _io_convert_ir_to_polynomials(
    origring::AbstractAlgebra.Generic.MPolyRing{T},
    gbexps::Vector{Vector{M}},
    gbcoeffs::Vector{Vector{C}},
    ord::Ord,
    input_ordering_matches_output
) where {T, M <: Monom, C <: Coeff, Ord <: AbstractMonomialOrdering}
    nv       = AbstractAlgebra.nvars(origring)
    ground   = AbstractAlgebra.base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    tmp      = Vector{aa_exponent_type}(undef, nv)
    @inbounds for i in 1:length(gbexps)
        cfs  = map(ground, gbcoeffs[i])
        exps = Vector{Vector{Int}}(undef, length(gbcoeffs[i]))
        for jt in 1:length(gbcoeffs[i])
            exps[jt] = gbexps[i][jt]
        end
        exported[i] = create_aa_polynomial(origring, cfs, exps)
    end
    exported
end
