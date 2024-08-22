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

function modular_reduce_mod_p!(
    ring::PolyRing,
    coeffs_zz::Vector{Vector{T1}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    coeffs_ff = [Vector{CoeffModular}(undef, length(c)) for c in coeffs_zz]
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

function modular_reduce_mod_p!(
    ring::PolyRing,
    basis::Basis{T1},
    prime::T2;
    deepcopy=true
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    ring_ff, coeffs_ff = modular_reduce_mod_p!(ring, basis.coeffs, prime)
    new_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(basis, coeffs_ff)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_ff)
    end
    ring_ff, new_basis
end

function modular_reduce_mod_p_in_batch!(
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

function modular_crt_full_changematrix!(
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
    end

    @inbounds for i in 1:length(state.changematrix_coeffs_zz)
        modulo = BigInt()
        crt_vec_full!(
            state.changematrix_coeffs_zz[i],
            modulo,
            [state.changematrix_coeffs_ff_all[ell][i] for ell in 1:length(lucky.used_primes)],
            map(eltype(eltype(eltype(state.changematrix_coeffs_ff_all[1]))), lucky.used_primes),
            [falses(length(state.changematrix_coeffs_ff_all[1][i][j])) for j in 1:length(state.changematrix_coeffs_ff_all[1][i])]
        )
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
    @inbounds for i in 1:length(state.changematrix_coeffs_zz)
        flag = ratrec_vec_full!(
            state.changematrix_coeffs_qq[i],
            state.changematrix_coeffs_zz[i],
            lucky.modulo,
            [falses(length(state.changematrix_coeffs_qq[i][j])) for j in 1:length(state.changematrix_coeffs_qq[i])]
        )
        !flag && return flag
    end
    true
end

function partial_rational_reconstruct!(
    state::ModularState,
    lucky::LuckyPrimes,
    witness_set::Vector{Tuple{Int, Int}}
)
    ratrec_vec_partial!(
        state.gb_coeffs_qq,
        state.gb_coeffs_zz,
        lucky.modulo,
        witness_set,
        state.ratrec_mask
    )
end

### 
# Checking correctness in modular computation

# Checks if the basis is reconstructed correctly.
# There are 3 levels of checks:
#   - heuristic check (dedicated to discard obviously bad cases),
#   - randomized check (checks correctness modulo a prime),
#   - certification (checks correctness directly over the rationals).
#
# Usually, by default, only the first two are active, which gives the correct
# basis with a high probability
function correctness_check!(
    state,
    lucky,
    ring,
    basis_qq,
    basis_zz,
    basis_ff,
    hashtable,
    params
)
    # First we check the size of the coefficients with a heuristic
    if params.heuristic_check
        if !modular_lift_heuristic_check(state.gb_coeffs_qq, lucky.modulo)
            @log :misc "Heuristic check failed."
            return false
        end
        @log :misc "Heuristic check passed!"
    end

    # Then check that a basis is also a basis modulo a prime
    if params.randomized_check
        if !randomized_correctness_check!(
            state,
            ring,
            basis_zz,
            basis_ff,
            lucky,
            hashtable,
            params
        )
            @log :misc "Randomized check failed."
            return false
        end
        @log :misc "Randomized check passed!"
    end
    if params.certify_check
        return certify_correctness_check!(
            state,
            ring,
            basis_qq,
            basis_ff,
            hashtable,
            params
        )
    end
    true
end

# Heuristic bound on the size of coefficients of the basis.
threshold_in_heuristic_check(num::BigInt, den::BigInt, modsz::Int) =
    1.10 * (Base.GMP.MPZ.sizeinbase(num, 2) + Base.GMP.MPZ.sizeinbase(den, 2)) >= modsz

# Checks that 
#   ln(num) + ln(den) < C ln(modulo)
# for all coefficients of form num/den
function modular_lift_heuristic_check(
    table_qq::Vector{Vector{T}},
    modulo::BigInt
) where {T <: CoeffQQ}
    modsz = Base.GMP.MPZ.sizeinbase(modulo, 2)
    @inbounds for i in 1:length(table_qq)
        for j in 1:length(table_qq[i])
            n = numerator(table_qq[i][j])
            d = denominator(table_qq[i][j])
            threshold_in_heuristic_check(n, d, modsz) && return false
        end
    end
    true
end

function modular_lift_heuristic_check_partial(
    table_qq::Vector{Vector{T}},
    modulo::BigInt,
    witness_set::Vector{Tuple{Int, Int}}
) where {T <: CoeffQQ}
    modsz = Base.GMP.MPZ.sizeinbase(modulo, 2)
    @inbounds for k in 1:length(witness_set)
        i, j = witness_set[k]
        n = numerator(table_qq[i][j])
        d = denominator(table_qq[i][j])
        threshold_in_heuristic_check(n, d, modsz) && return false
    end
    true
end

function randomized_correctness_check!(
    state,
    ring,
    input_zz,
    gb_ff,
    lucky,
    hashtable,
    params
)
    # !!! note that this function may modify the given hashtable!
    prime = primes_next_aux_prime!(lucky)
    ring_ff, input_ff = modular_reduce_mod_p!(ring, input_zz, prime, deepcopy=true)
    # TODO: do we really need to re-scale things to be fraction-free?
    gb_coeffs_zz = _clear_denominators!(state.gb_coeffs_qq)
    gb_zz = basis_deep_copy_with_new_coeffs(gb_ff, gb_coeffs_zz)
    ring_ff, gb_ff = modular_reduce_mod_p!(ring, gb_zz, prime, deepcopy=false)
    arithmetic = select_arithmetic(CoeffModular, prime, :auto, false)
    basis_make_monic!(gb_ff, arithmetic, params.changematrix)
    # Check that some polynomial is not reduced to zero
    f4_normalform!(ring_ff, gb_ff, input_ff, hashtable, arithmetic)
    for i in 1:(input_ff.nprocessed)
        if !io_iszero_coeffs(input_ff.coeffs[i])
            return false
        end
    end
    # Check that the basis is a groebner basis
    pairset = pairset_initialize(UInt64)
    if !f4_isgroebner!(ring_ff, gb_ff, pairset, hashtable, arithmetic)
        return false
    end
    true
end

function certify_correctness_check!(state, ring, input_qq, gb_ff, hashtable, params)
    gb_qq = basis_deep_copy_with_new_coeffs(gb_ff, state.gb_coeffs_qq)
    ring_qq = PolyRing(ring.nvars, ring.ord, 0)
    input_qq = basis_deepcopy(input_qq)
    # Check that some polynomial is not reduced to zero modulo a prime
    f4_normalform!(ring_qq, gb_qq, input_qq, hashtable, params.arithmetic)
    for i in 1:(input_qq.nprocessed)
        if !io_iszero_coeffs(input_qq.coeffs[i])
            return false
        end
    end
    # Check that the basis is a groebner basis modulo a prime
    pairset = pairset_initialize(UInt64)
    if !f4_isgroebner!(ring_qq, gb_qq, pairset, hashtable, params.arithmetic)
        return false
    end
    true
end

# TODO :)
function majority_vote!(state, basis_ff, params)
    true
end
