
#------------------------------------------------------------------------------

mutable struct PolyRing
    #= Ring information =#
    # number of variables
    nvars::Int
    # raw length of exponent vector
    explen::Int
    # ring monomial ordering,
    # possible are :lex and :degrevlex
    ord::Symbol
    # characteristic of coefficient field
    ch::UInt64
end

#------------------------------------------------------------------------------

# Finite field AbstractAlgebra.MPoly conversion specialization
#
# converts MPoly representation to internal polynomial representation
# extracting base ring, exponents, and coefficients
function convert_to_internal(orig_polys::Vector{MPoly{GFElem{Int64}}})
    isempty(orig_polys) && error("Empty groebner input")

    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            # TODO: ask
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = map(UInt64 ∘ data, coefficients(poly))
    end

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord
    ch     = characteristic(R)

    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    @assert 0 < ch < 2^32
    # @assert ord == :degrevlex
    # @assert nvars > 1 && nvars + 1 == explen

    return ring, exps, coeffs
end

# Rational field AbstractAlgebra.MPoly conversion specialization
#
# converts MPoly representation to internal polynomial representation
# extracting base ring, exponents, and coefficients
function convert_to_internal(orig_polys::Vector{MPoly{Rational{BigInt}}})
    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            # TODO: ask
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = copy(poly.coeffs)
    end

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord
    ch     = 0

    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    # @assert ord == :degrevlex
    # @assert nvars > 1 && nvars + 1 == explen

    return ring, exps, coeffs
end

#------------------------------------------------------------------------------

# Finite field AbstractAlgebra.MPoly export specialization
#
function export_basis(
            origring::MPolyRing{GFElem{Int64}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(gbexps[i]...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

# Rational field AbstractAlgebra.MPoly export specialization
#
function export_basis(
            origring::MPolyRing{Rational{BigInt}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{Rational{BigInt}}})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(gbexps[i]...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------
