# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Rational reconstruction

# Returns the largest integer N (possibly off by 1) such that 2 N^2 < m.
function ratrec_reconstruction_bound(m::BigInt)
    isqrt((m >> UInt64(1)) - 1)
end

function ratrec_nemo(
    a::Nemo.ZZRingElem,
    m::Nemo.ZZRingElem,
    N::Nemo.ZZRingElem,
    D::Nemo.ZZRingElem
)
    success, pq = Nemo.reconstruct(a, m, N, D)
    success, Rational{BigInt}(pq)
end

###
# Element-wise rational reconstruction

# table of rationals need not be initialized
function ratrec_vec_partial!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    witness_set::Vector{Tuple{Int, Int}},
    mask::Vector{BitVector}
)
    @invariant length(table_qq) == length(table_zz)

    bound = ratrec_reconstruction_bound(modulo)
    nemo_bound = Nemo.ZZRingElem(bound)
    nemo_modulo = Nemo.ZZRingElem(modulo)

    @inbounds for k in 1:length(witness_set)
        i, j = witness_set[k]
        @invariant 1 <= i <= length(table_zz) && 1 <= j <= length(table_zz[i])

        rem_nemo = Nemo.ZZRingElem(table_zz[i][j])

        success, pq = ratrec_nemo(rem_nemo, nemo_modulo, nemo_bound, nemo_bound)
        !success && return false

        table_qq[i][j] = pq

        mask[i][j] = true
    end

    true
end

function ratrec_vec_full!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    mask::Vector{BitVector}
)
    @invariant length(table_qq) == length(table_zz)
    @invariant modulo > 1

    bound = ratrec_reconstruction_bound(modulo)
    nemo_bound = Nemo.ZZRingElem(bound)
    nemo_modulo = Nemo.ZZRingElem(modulo)

    @inbounds for i in 1:length(table_zz)
        @invariant length(table_zz[i]) == length(table_qq[i])
        for j in 1:length(table_zz[i])
            if mask[i][j]
                continue
            end

            rem_nemo = Nemo.ZZRingElem(table_zz[i][j])
            @invariant 0 <= rem_nemo < modulo

            success, pq = ratrec_nemo(rem_nemo, nemo_modulo, nemo_bound, nemo_bound)
            !success && return false

            table_qq[i][j] = pq
        end
    end

    true
end
