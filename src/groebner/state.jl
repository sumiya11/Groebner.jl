# This file is a part of Groebner.jl. License is GNU GPL v2.

mutable struct ModularState{T1 <: CoeffZZ, T2 <: CoeffQQ, T3}
    gb_coeffs_zz::Vector{Vector{T1}}
    gb_coeffs_qq::Vector{Vector{T2}}
    gb_coeffs_ff_all::Vector{Vector{Vector{T3}}}

    crt_mask::Vector{BitVector}
    ratrec_mask::Vector{BitVector}

    changematrix_coeffs_ff_all::Vector{Vector{Vector{Vector{T3}}}}
    changematrix_coeffs_zz::Vector{Vector{Vector{T1}}}
    changematrix_coeffs_qq::Vector{Vector{Vector{T2}}}

    function ModularState{T1, T2, T3}(
        params::AlgorithmParameters
    ) where {T1 <: CoeffZZ, T2 <: CoeffQQ, T3 <: CoeffZp}
        new(
            Vector{Vector{T1}}(),
            Vector{Vector{T2}}(),
            Vector{Vector{Vector{T3}}}(),
            Vector{BitVector}(),
            Vector{BitVector}(),
            Vector{Vector{Vector{Vector{T3}}}}(),
            Vector{Vector{Vector{T1}}}(),
            Vector{Vector{Vector{T2}}}()
        )
    end
end

function modular_witness_set(gb_coeffs, params::AlgorithmParameters)
    witness_set = Vector{Tuple{Int, Int}}(undef, length(gb_coeffs))
    rng = params.rng
    k = 1
    @inbounds for i in 1:length(gb_coeffs)
        l = length(gb_coeffs[i])
        isone(l) && continue
        nl = max(floor(Int, log2(l)) - 1, 1)
        while k + nl + 1 > length(witness_set)
            resize!(witness_set, 2 * length(witness_set))
        end
        for j in 1:nl
            witness_set[k] = (i, rand(rng, 2:l))
            k += 1
        end
        witness_set[k] = (i, l)
        k += 1
    end
    resize!(witness_set, k - 1)
    unique!(witness_set)
    witness_set
end

function common_denominator!(den::BigInt, coeffs::Vector{T}) where {T <: CoeffQQ}
    Base.GMP.MPZ.set_si!(den, 1)
    @inbounds for i in 1:length(coeffs)
        c = coeffs[i]
        Base.GMP.MPZ.lcm!(den, denominator(c))
    end
    den
end

function common_denominator(coeffs::Vector{T}) where {T <: CoeffQQ}
    common_denominator!(BigInt(), coeffs)
end

function _clear_denominators!(coeffs_qq::Vector{Vector{T}}) where {T <: CoeffQQ}
    coeffs_zz = [[BigInt(0) for _ in 1:length(c)] for c in coeffs_qq]
    den, buf = BigInt(), BigInt()
    @inbounds for i in 1:length(coeffs_qq)
        @invariant length(coeffs_zz[i]) == length(coeffs_qq[i])
        den = common_denominator!(den, coeffs_qq[i])
        for j in 1:length(coeffs_qq[i])
            num = numerator(coeffs_qq[i][j])
            Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs_qq[i][j]))
            Base.GMP.MPZ.mul!(coeffs_zz[i][j], num, buf)
        end
    end
    coeffs_zz
end

function clear_denominators!(basis::Basis{T}; deepcopy=false) where {T <: CoeffQQ}
    coeffs_zz = _clear_denominators!(basis.coeffs)
    if deepcopy
        basis_deep_copy_with_new_coeffs(basis, coeffs_zz)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_zz)
    end
end

function bigint_mod_p!(buf::BigInt, x::BigInt, prime::Unsigned, prime_big::BigInt)
    if Base.GMP.MPZ.cmp_ui(x, 0) < 0
        Base.GMP.MPZ.fdiv_q!(buf, x, prime_big)
        Base.GMP.MPZ.neg!(buf)
        Base.GMP.MPZ.mul_ui!(buf, prime)
        Base.GMP.MPZ.add!(x, buf)
    end
    Base.GMP.MPZ.tdiv_r!(buf, x, prime_big)
    buf
end

function reduce_modulo_p!(
    ring,
    coeffs_zz::Vector{Vector{T1}},
    coeffs_ff::Vector{Vector{T2}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    p   = BigInt()
    buf = BigInt()
    c   = BigInt()
    Base.GMP.MPZ.set_ui!(p, prime)

    @inbounds for i in 1:length(coeffs_zz)
        cfs_zz_i = coeffs_zz[i]
        for j in 1:length(cfs_zz_i)
            Base.GMP.MPZ.set!(c, cfs_zz_i[j])
            bigint_mod_p!(buf, c, prime, p)
            coeffs_ff[i][j] = CoeffModular(buf)
        end
    end
    ring_ff = PolyRing(ring.nvars, ring.ord, UInt(prime))
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    ring::PolyRing,
    coeffs_zz::Vector{Vector{T1}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    coeffs_ff = [Vector{CoeffModular}(undef, length(c)) for c in coeffs_zz]
    ring_ff, coeffs_ff = reduce_modulo_p!(ring, coeffs_zz, coeffs_ff, prime)
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    ring::PolyRing,
    basis::Basis{T1},
    prime::T2;
    deepcopy=true
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    ring_ff, coeffs_ff = reduce_modulo_p!(ring, basis.coeffs, prime)
    new_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(basis, coeffs_ff)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_ff)
    end
    ring_ff, new_basis
end

function reduce_modulo_p_in_batch!(
    ring::PolyRing,
    basis::Basis{C},
    prime_xn::NTuple{N, T}
) where {C, N, T}
    coeffs_zz = basis.coeffs
    coeffs_ff_xn = [Vector{CompositeNumber{N, T}}(undef, length(c)) for c in coeffs_zz]

    p = BigInt()
    buf = BigInt()
    xn = map(_ -> BigInt(0), 1:N)
    c = BigInt()
    prime_big_xn = map(BigInt, prime_xn)

    @inbounds for i in 1:length(coeffs_zz)
        cfs_zz_i = coeffs_zz[i]
        for j in 1:length(cfs_zz_i)
            for k in 1:N
                Base.GMP.MPZ.set!(xn[k], cfs_zz_i[j])
            end
            data = ntuple(
                k -> T(bigint_mod_p!(buf, xn[k], UInt(prime_xn[k]), prime_big_xn[k])),
                N
            )
            coeffs_ff_xn[i][j] = CompositeNumber{N, T}(data)
        end
    end
    ring_ff_4x = PolyRing(ring.nvars, ring.ord, CompositeNumber{N, T}(prime_xn))
    basis_ff_4x = basis_deep_copy_with_new_coeffs(basis, coeffs_ff_xn)

    ring_ff_4x, basis_ff_4x
end

function resize_state_if_needed!(
    state::ModularState,
    gb_coeffs::Vector{Vector{T}}
) where {T}
    resize!(state.gb_coeffs_zz, length(gb_coeffs))
    resize!(state.gb_coeffs_qq, length(gb_coeffs))
    resize!(state.crt_mask, length(gb_coeffs))
    resize!(state.ratrec_mask, length(gb_coeffs))

    @inbounds for i in 1:length(gb_coeffs)
        state.gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.gb_coeffs_qq[i] = [Rational{BigInt}(1) for _ in 1:length(gb_coeffs[i])]
        state.crt_mask[i] = falses(length(gb_coeffs[i]))
        state.ratrec_mask[i] = falses(length(gb_coeffs[i]))
    end

    nothing
end

###
# CRT

function modular_crt_full!(state::ModularState, lucky::LuckyPrimes)
    if isempty(state.gb_coeffs_zz)
        resize_state_if_needed!(state, state.gb_coeffs_ff_all[1])
    end
    crt_vec_full!(
        state.gb_coeffs_zz,
        lucky.modulo,
        state.gb_coeffs_ff_all,
        map(eltype(eltype(state.gb_coeffs_ff_all[1])), lucky.used_primes),
        state.crt_mask
    )
    @invariant lucky.modulo == prod(BigInt, lucky.used_primes)
    nothing
end

function full_simultaneous_crt_reconstruct_changematrix!(
    state::ModularState,
    lucky::LuckyPrimes
)
    if isempty(state.changematrix_coeffs_zz)
        coeffs_ff = state.changematrix_coeffs_ff_all[1]

        resize!(state.changematrix_coeffs_zz, length(coeffs_ff))
        resize!(state.changematrix_coeffs_qq, length(coeffs_ff))

        @inbounds for i in 1:length(coeffs_ff)
            state.changematrix_coeffs_zz[i] =
                Vector{Vector{BigInt}}(undef, length(coeffs_ff[i]))
            state.changematrix_coeffs_qq[i] =
                Vector{Vector{Rational{BigInt}}}(undef, length(coeffs_ff[i]))
            for j in 1:length(coeffs_ff[i])
                state.changematrix_coeffs_zz[i][j] =
                    [BigInt(0) for _ in 1:length(coeffs_ff[i][j])]
                state.changematrix_coeffs_qq[i][j] =
                    [Rational{BigInt}(1) for _ in 1:length(coeffs_ff[i][j])]
            end
        end

        changematrix_coeffs_zz = state.changematrix_coeffs_zz
        @inbounds for i in 1:length(coeffs_ff)
            for j in 1:length(coeffs_ff[i])
                for k in 1:length(coeffs_ff[i][j])
                    Base.GMP.MPZ.set_ui!(
                        changematrix_coeffs_zz[i][j][k],
                        coeffs_ff[i][j][k]
                    )
                end
            end
        end
        return nothing
    end

    changematrix_coeffs_zz = state.changematrix_coeffs_zz

    # Takes the lock..
    @invariant length(state.changematrix_coeffs_ff_all) == length(lucky.used_primes)

    n = length(lucky.used_primes)
    rems = Vector{UInt64}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.used_primes
    crt_precompute!(M, n1, n2, mults, moduli)

    @inbounds for i in 1:length(changematrix_coeffs_zz)
        for j in 1:length(changematrix_coeffs_zz[i])
            for k in 1:length(changematrix_coeffs_zz[i][j])
                for ell in 1:length(lucky.used_primes)
                    rems[ell] = state.changematrix_coeffs_ff_all[ell][i][j][k] % UInt64
                end
                crt!(M, buf, n1, n2, rems, mults)
                Base.GMP.MPZ.set!(changematrix_coeffs_zz[i][j][k], buf)
            end
        end
    end

    nothing
end

function modular_crt_partial!(
    state::ModularState,
    lucky::LuckyPrimes,
    witness_set::Vector{Tuple{Int, Int}}
)
    crt_vec_partial!(
        state.gb_coeffs_zz,
        lucky.modulo,
        state.gb_coeffs_ff_all,
        map(eltype(eltype(state.gb_coeffs_ff_all[1])), lucky.used_primes),
        witness_set,
        state.crt_mask
    )
    nothing
end

function modular_ratrec_vec_full!(state::ModularState, lucky::LuckyPrimes)
    @invariant lucky.modulo == prod(BigInt, lucky.used_primes)
    ratrec_vec_full!(
        state.gb_coeffs_qq,
        state.gb_coeffs_zz,
        lucky.modulo,
        state.ratrec_mask
    )
end

function full_rational_reconstruct_changematrix!(state::ModularState, lucky::LuckyPrimes)
    modulo = lucky.modulo
    @invariant modulo == prod(BigInt, lucky.used_primes)

    bnd = ratrec_reconstruction_bound(modulo)

    changematrix_coeffs_zz = state.changematrix_coeffs_zz
    changematrix_coeffs_qq = state.changematrix_coeffs_qq
    # ratrec_mask = state.ratrec_mask

    @invariant length(changematrix_coeffs_zz) == length(changematrix_coeffs_qq)
    nemo_modulo = Nemo.ZZRingElem(modulo)
    nemo_bnd = Nemo.ZZRingElem(bnd)

    @inbounds for i in 1:length(changematrix_coeffs_zz)
        @invariant length(changematrix_coeffs_zz[i]) == length(changematrix_coeffs_qq[i])
        for j in 1:length(changematrix_coeffs_zz[i])
            for k in 1:length(changematrix_coeffs_zz[i][j])
                cz = changematrix_coeffs_zz[i][j][k]
                nemo_rem = Nemo.ZZRingElem(cz)
                success, pq = ratrec_nemo_2(nemo_rem, nemo_modulo, nemo_bnd, nemo_bnd)
                !success && return false

                changematrix_coeffs_qq[i][j][k] = pq
            end
        end
    end

    true
end

function partial_rational_reconstruct!(
    state::ModularState,
    lucky::LuckyPrimes,
    witness_set::Vector{Tuple{Int, Int}}
)
    ratrec_vec_partial!(state.gb_coeffs_qq, state.gb_coeffs_zz, lucky.modulo, witness_set, state.ratrec_mask)
end
