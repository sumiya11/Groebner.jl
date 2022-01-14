
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
end

#------------------------------------------------------------------------------

"""
    Converts input polynomials to internal representation used by the algorithm.
    Extracts base ring information, exponents, and coefficients.

    Currently, there are specializations for:
        `AbstractAlgebra.MPoly` polynomials
"""
function convert_to_internal(orig_polys::Vector{T}) where {T}
    isempty(orig_polys) && error("Empty input")

    convert_to_internal(orig_polys)
end

"""
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(orig_polys::Vector{MPoly{T}}) where {T}
    isempty(orig_polys) && error("Empty input")

    ord = ordering(parent(first(orig_polys)))
    convert_to_internal(orig_polys, Val(ord))
end

"""
    Finite field :degrevlex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{MPoly{GFElem{Int64}}},
        ::Val{:degrevlex})

    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord
    ch     = characteristic(R)

    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    @assert 0 < ch < 2^32
    @assert ord == :degrevlex
    @assert nvars > 1 && nvars + 1 == explen

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = map(UInt64 ∘ data, coefficients(poly))
    end

    return ring, exps, coeffs
end

"""
    Finite field :lex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{MPoly{GFElem{Int64}}},
        ::Val{:lex})

    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    R = parent(first(orig_polys))
    explen = R.N + 1
    nvars  = R.num_vars
    ord    = R.ord
    ch     = characteristic(R)
    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    @assert 0 < ch < 2^32
    @assert ord == :lex
    @assert nvars > 1 && nvars + 1 == explen

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, explen)
            exps[i][j][end] = zero(exps[i][j][end])
            for ei in 1:nvars
                exps[i][j][ei] = poly.exps[nvars - ei + 1, j]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
        coeffs[i] = map(UInt64 ∘ data, coefficients(poly))
    end

    return ring, exps, coeffs
end

"""
    Rational field :degrevlex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{MPoly{Rational{BigInt}}},
        ::Val{:degrevlex})

    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord
    ch     = 0
    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    @assert ord == :degrevlex
    @assert nvars > 1 && nvars + 1 == explen

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = copy(poly.coeffs)
    end

    return ring, exps, coeffs
end


"""
    Rational field :lex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function convert_to_internal(
        orig_polys::Vector{MPoly{Rational{BigInt}}},
        ::Val{:lex})

    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord
    ch     = 0
    ring = PolyRing(nvars, explen, ord, UInt64(ch))

    @assert ord == :degrevlex
    @assert nvars > 1 && nvars + 1 == explen

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, explen)
            exps[i][j][end] = zero(exps[i][j][end])
            for ei in 1:nvars
                exps[i][j][ei] = poly.exps[nvars - ei + 1, j]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
        coeffs[i] = copy(poly.coeffs)
    end

    return ring, exps, coeffs
end

#------------------------------------------------------------------------------

"""
    Converts internal polynomials for export as elements of `origring`.

    Currently, there are specializations for:
        `AbstractAlgebra.MPolyRing`
"""
function export_basis(
            origring::T,
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}}) where {T}

    export_basis(origring, gbexps, gbcoeffs)
end

"""
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            origring::MPolyRing{T},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}}) where {T}

    ord = ordering(origring)
    export_basis(origring, gbexps, gbcoeffs, Val(ord))
end

"""
    Finite field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            origring::MPolyRing{GFElem{Int64}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}},
            ::Val{:degrevlex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(gbexps[i]...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

"""
    Finite field :lex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            origring::MPolyRing{GFElem{Int64}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}},
            ::Val{:lex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> x[end-1:-1:1], gbexps[i])...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

"""
    Rational field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            origring::MPolyRing{Rational{BigInt}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{Rational{BigInt}}},
            ::Val{:degrevlex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(gbexps[i]...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

"""
    Rational field :lex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            origring::MPolyRing{Rational{BigInt}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{Rational{BigInt}}},
            ::Val{:lex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> x[end-1:-1:1], gbexps[i])...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------
