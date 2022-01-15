
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
        # convert_to_internal(orig_polys, ordering, Val(:noparent))
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

function extract_ring(
        orig_polys::Vector{<:AbstractPolynomial{T}}) where {T <: Union{Integer, Rational}}

    f = first(orig_polys)

    nv = maximum(map(FastGroebner.MultivariatePolynomials.nvariables, orig_polys))
    explen = nv + 1
    ord    = :deglex
    ch     = 0

    @assert nv + 1 == explen
    @assert ord in (:lex, :degrevlex, :deglex)

    PolyRing(nv, explen, ord, UInt64(ch), :multivariate)
end


function extract_ring(
        orig_polys::Vector{<:AbstractPolynomial{T}}) where {T <: AbstractAlgebra.Generic.FieldElem}

    f = first(orig_polys)

    nv = maximum(map(FastGroebner.MultivariatePolynomials.nvariables, orig_polys))
    explen = nv + 1
    ord    = :deglex
    ch     = characteristic(parent(first(coefficients(f))))

    @assert nv + 1 == explen
    @assert ord in (:lex, :degrevlex, :deglex)

    PolyRing(nv, explen, ord, UInt64(ch), :multivariate)
end

function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{<:FastGroebner.AbstractPolynomial{T}}) where {T <: AbstractAlgebra.Generic.FieldElem}

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    explen = ring.explen
    nvars  = ring.nvars

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for (j, t) in zip(1:length(poly), MultivariatePolynomials.terms(poly))
            exps[i][j] = Vector{UInt16}(undef, explen)
            exps[i][j][end] = zero(exps[i][j][end])
            et = MultivariatePolynomials.exponents(t)
            for ei in 1:nvars
                exps[i][j][ei] = et[ei]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
        coeffs[i] = map(UInt64 ∘ data, MultivariatePolynomials.coefficients(poly))
    end

    exps, coeffs
end

function exponents_wrt_vars(t, var2idx)
    exp = zeros(UInt16, length(var2idx))
    for (v, p) in FastGroebner.MultivariatePolynomials.powers(t)
        exp[var2idx[v]] = p
    end
    exp
end

function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{<:FastGroebner.AbstractPolynomial{T}}) where {T <: Union{Integer, Rational}}

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    explen = ring.explen
    nvars  = ring.nvars
    vars   = union!(map(FastGroebner.MultivariatePolynomials.variables, orig_polys)...)
    sort!(vars, rev=true)
    var2idx = Dict(vars[i] => i for i in 1:length(vars))

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for (j, t) in enumerate(FastGroebner.MultivariatePolynomials.monomials(poly))
            exps[i][j] = Vector{UInt16}(undef, explen)
            exps[i][j][end] = zero(exps[i][j][end])
            et = exponents_wrt_vars(t, var2idx)
            for ei in 1:nvars
                exps[i][j][ei] = et[ei]
                exps[i][j][end] += exps[i][j][ei]
            end
        end
        coeffs[i] = map(Rational, FastGroebner.MultivariatePolynomials.coefficients(poly))
    end

    exps, coeffs
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

"""
    Finite field :lex and :deglex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{MPoly{GFElem{Int64}}},
            ::Ord) where {Ord<:Union{Val{:lex}, Val{:deglex}}}

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    explen = ring.explen
    nvars  = ring.nvars

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

    exps, coeffs
end

"""
    Finite field :degrevlex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{MPoly{GFElem{Int64}}},
            ::Val{:degrevlex})

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{UInt64}}(undef, npolys)

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = map(UInt64 ∘ data, coefficients(poly))
    end

    exps, coeffs
end

"""
    Rational field :degrevlex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{MPoly{Rational{BigInt}}},
            ::Val{:degrevlex})

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = copy(poly.coeffs)
    end

    exps, coeffs
end

"""
    Rational field :lex and :deglex
    `AbstractAlgebra.MPoly` conversion specialization
"""
function extract_polys(
            ring::PolyRing,
            orig_polys::Vector{MPoly{Rational{BigInt}}},
            ::Ord) where {Ord<:Union{Val{:lex}, Val{:deglex}}}

    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Rational{BigInt}}}(undef, npolys)

    explen = ring.explen
    nvars  = ring.nvars

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

    exps, coeffs
end

function extract_ring(R::MPolyRing{T}) where {T}
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

    PolyRing(nv, explen, ord, UInt64(ch), :abstract)
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
    exps, cfs = extract_polys(ring, orig_polys, Val(ring.ord))
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
    elseif ring.origring == :dynamic
        convert_to_output(ring, origpolys, gbexps, gbcoeffs)
    elseif ring.origring == :hasparent
        convert_to_output(ring, parent(first(origpolys)), gbexps, gbcoeffs)
    else
        error("Sorry, unknown polynomial ring.")
    end
end

#------------------------------------------------------------------------------

"""
    `multivariate` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {P<:AbstractPolynomial{J}, I<:Rational} where {J}

    ground = I
    origvars = MultivariatePolynomials.variables(first(origpolys))
    exported = Vector{P}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        expvectors = [gbexps[i][j][1:end-1] for j in 1:length(gbexps[i])]
        expvars = map(ev -> prod(origvars .^ ev), expvectors)
        exported[i] = P(cfs, expvars)
    end
    exported
end


"""
    `multivariate` conversion specialization
"""
function convert_to_output(
            ring::PolyRing,
            origpolys::Vector{P},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {P<:AbstractPolynomial, I<:AbstractAlgebra.Generic.FieldElem}

    ground = parent(first(coefficients(first(origpolys))))
    origvars = variables(first(origpolys))
    exported = Vector{P}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        expvectors = [gbexps[i][j][1:end-1] for j in 1:length(gbexps[i])]
        expvars = map(ev -> prod(origvars .^ ev), expvectors)
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

    convert_to_output(origring, gbexps, gbcoeffs)
end

function convert_to_output(
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
function convert_to_output(
            ring::PolyRing,
            origring::MPolyRing{T},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{I}}) where {T, I}

    ord = ordering(origring)
    convert_to_output(origring, gbexps, gbcoeffs, Val(ord))
end

"""
    Finite field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
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
function convert_to_output(
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
    Finite field :deglex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{GFElem{Int64}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{UInt64}},
            ::Val{:deglex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> [x[end-1:-1:1]..., x[end]], gbexps[i])...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

"""
    Rational field :degrevlex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
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
function convert_to_output(
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

"""
    Rational field :deglex
    `AbstractAlgebra.MPolyRing` conversion specialization
"""
function convert_to_output(
            origring::MPolyRing{Rational{BigInt}},
            gbexps::Vector{Vector{Vector{UInt16}}},
            gbcoeffs::Vector{Vector{Rational{BigInt}}},
            ::Val{:deglex})

    ground   = base_ring(origring)
    exported = Vector{elem_type(origring)}(undef, length(gbexps))
    for i in 1:length(gbexps)
        cfs    = map(ground, gbcoeffs[i])
        exps   = UInt64.(hcat(map(x -> [x[end-1:-1:1]..., x[end]], gbexps[i])...))
        exported[i] = MPoly{elem_type(ground)}(origring, cfs, exps)
    end
    exported
end

#------------------------------------------------------------------------------
