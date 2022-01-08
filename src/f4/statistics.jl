
#=
    The file for keeping track of intermediate results.
    It mirrors some f4 functions to collect statistics
    while computing the basis.
=#

#------------------------------------------------------------------------------

# Contains information about one iteration
mutable struct IterInfo
    basissize::Int
    pairsetsize::Int
    pairsselected::Int
    matrixsize::Tuple{Int, Int}
    matrixdensity::Float64
    newpivots::Int
end

# Contains information about input generators
struct BasisInfo
    ngens::Int
    maxtotaldeg::Int

    function BasisInfo()
        new(0, 0)
    end
    function BasisInfo(ngens::Int, maxtotaldeg::Int)
        new(ngens, maxtotaldeg)
    end
end

# Stores all collected statistics for one groebner computation
mutable struct Statistics
    inputs::BasisInfo
    output::BasisInfo

    data::Vector{IterInfo}

    function Statistics()
        new(BasisInfo(), BasisInfo(), Vector{IterInfo}())
    end
end

#------------------------------------------------------------------------------

function collect_stats(basis::Basis)
    ngens  = basis.ntotal
    # maxdeg = max(last.(basis.exponents))
    maxdeg = 0
    BasisInfo(ngens, maxdeg)
end

#------------------------------------------------------------------------------

function groebner(
            polys::Vector{MPoly{GFElem{Int}}},
            stats::Statistics;
            rng::Rng=Random.MersenneTwister(42)
            ) where {Rng<:Random.AbstractRNG}

    length(polys) > 0 || error("Empty input")
    R = parent(first(polys))
    ordering(R) in (:degrevlex, ) || error("Only :degrevlex is supported")

    ring, exps, coeffs = convert_to_internal(polys)

    gb, ht = f4(ring, exps, coeffs, rng, stats)

    export_basis(parent(first(polys)), gb, ht)
end

function f4(ring::PolyRing,
            exponents::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            rng::Rng,
            stats::Statistics;
            reduced::Bool=true,
            tablesize::Int=2^16
            ) where {Rng<:Random.AbstractRNG}

    gb, ht = f4(ring, exponents, coeffs, rng; reduced=reduced, tablesize=tablesize)

    stats.output = collect_stats(gb)

    gb, ht
end

function initialize_structures(
            ring::PolyRing,
            exponents::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            rng::Random.AbstractRNG,
            tablesize::Int,
            stats::Statistics)

    basis, ht = initialize_structures(ring, exponents, coeffs, rng, tablesize)

    stats.inputs = collect_stats(basis)

    basis, ht
end
