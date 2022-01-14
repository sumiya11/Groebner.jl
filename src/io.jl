
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
    #    :dynamic for DynamicPolynomials,
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
"""
function convert_to_internal(orig_polys::Vector{T}) where {T}
    isempty(orig_polys) && error("Empty input")

    if hasmethod(parent, Tuple{typeof(first(orig_polys))})
        convert_to_internal(orig_polys, Val(:hasparent))
    else
        error("Sorry, we don't work with this type of polynomials yet. Feel free to open a github issue")
        # convert_to_internal(orig_polys, Val(:noparent))
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
    @assert ord in (:lex, :degrevlex)
    @assert 0 <= ch < 2^32

    PolyRing(nv, explen, ord, UInt64(ch), :hasparent)
end

function extract_polys(ring::PolyRing, orig_polys::Vector{T}) where {T}
    if ring.ch > 0
        extract_polys_ff(ring, orig_polys)
    else
        extract_polys_qq(ring, orig_polys)
    end
end

function extract_polys_ff(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    f = first(orig_polys)
    @assert hasmethod(exponent_vector, Tuple{typeof(f)})
    @assert hasmethod(coefficients, Tuple{typeof(f)})

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, ring.explen)
            exps[i][j][1:ring.nvars] .= exponent_vector(poly, j)
            exps[i][j][end] = sum(exps[i][j][k] for k in 1:ring.nvars)
        end
        coeffs[i] = map(UInt64 ∘ data, coefficients(poly))
    end

    exps, coeffs
end

function extract_polys_qq(ring::PolyRing, orig_polys::Vector{T}) where {T}
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    f = first(orig_polys)
    @assert hasmethod(exponent_vector, Tuple{typeof(f), Int})
    @assert hasmethod(coefficients, Tuple{typeof(f)})
    @assert hasmethod(Rational, Tuple{typeof(leading_coefficient(f))})

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = Vector{UInt16}(undef, ring.explen)
            exps[i][j][1:ring.nvars] .= exponent_vector(poly, j)
            exps[i][j][end] = sum(exps[i][j][k] for k in 1:ring.nvars)
        end
        coeffs[i] = map(Rational, coefficients(poly))
    end

    exps, coeffs
end

function convert_to_internal(
            orig_polys::Vector{T},
            ::Val{:hasparent}) where {T}

    R = parent(first(orig_polys))
    ring = extract_ring(R)
    exps, cfs = extract_polys(ring, orig_polys)

    ring, exps, cfs
end

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

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

    ring = PolyRing(nvars, explen, ord, UInt64(ch), :abstract)

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
    ring = PolyRing(nvars, explen, ord, UInt64(ch), :abstract)

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
    ring = PolyRing(nvars, explen, ord, UInt64(ch), :abstract)

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
    explen = R.N + 1
    nvars  = R.num_vars
    ord    = R.ord
    ch     = 0
    ring = PolyRing(nvars, explen, ord, UInt64(ch), :abstract)

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
        coeffs[i] = copy(poly.coeffs)
    end

    return ring, exps, coeffs
end

#------------------------------------------------------------------------------

"""
    Converts internal polynomials for export as elements of `origring`.

    This is the most general implementation.
    It happened to work for polynomials from
        . `Nemo`

    Currently, there are more efficient specializations for:
        . `AbstractAlgebra.MPoly`
"""
function export_basis(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {P, I}

    if ring.origring == :abstract
        export_basis(ring, parent(first(origpolys)), gbexps, gbcoeffs)
    elseif ring.origring == :dynamic
        # export_basis(ring, parent(first(origpolys)), gbexps, gbcoeffs)
    elseif ring.origring == :hasparent
        export_basis(ring, parent(first(origpolys)), gbexps, gbcoeffs)
    else
        error("Sorry, unknown polynomial ring.")
    end
end

#------------------------------------------------------------------------------

"""
    `hasparent` conversion specialization
"""
function export_basis(
            ring::PolyRing,
            origring::M,
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {M, T, I}

    export_basis(origring, gbexps, gbcoeffs)
end

function export_basis(
            origring::M,
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where{M, I}

    @assert hasmethod(base_ring, Tuple{typeof(origring)})

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

"""
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function export_basis(
            ring::PolyRing,
            origring::MPolyRing{T},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {T, I}

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
