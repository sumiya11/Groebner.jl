# This file is a part of Groebner.jl. License is GNU GPL v2.

# Acts as a buffer when doing arithmetic with big rationals
mutable struct CoefficientBuffer
    # Buffers for scaling coefficients
    scalebuf1::BigInt
    scalebuf2::BigInt
    # Buffers for reductions modulo a prime
    reducebuf1::BigInt
    reducebuf2::BigInt
    reducebuf3::BigInt
    # Buffers for coefficient reconstruction
    reconstructbuf1::BigInt
    reconstructbuf2::BigInt
    reconstructbuf3::BigInt
    reconstructbuf4::BigInt
    reconstructbuf5::BigInt
    reconstructbuf6::BigInt
    reconstructbuf7::BigInt
    reconstructbuf8::BigInt
    reconstructbuf9::BigInt
    reconstructbuf10::BigInt

    function CoefficientBuffer()
        new(
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt(),
            BigInt()
        )
    end
end

# The state of the modular GB 
mutable struct GroebnerState{T1 <: CoeffZZ, T2 <: CoeffQQ, T3}
    gb_coeffs_zz::Vector{Vector{T1}}
    gb_coeffs_qq::Vector{Vector{T2}}

    gb_coeffs_ff_all::Vector{Vector{Vector{T3}}}
    prev_index::Int

    selected_coeffs_zz::Vector{T1}
    selected_prev_coeffs_zz::Vector{T1}
    selected_coeffs_qq::Vector{T2}
    is_crt_reconstructed_mask::Vector{BitVector}
    is_rational_reconstructed_mask::Vector{BitVector}

    buffer::CoefficientBuffer

    changematrix_coeffs_ff_all::Vector{Vector{Vector{Vector{T3}}}}
    changematrix_coeffs_zz::Vector{Vector{Vector{T1}}}
    changematrix_coeffs_qq::Vector{Vector{Vector{T2}}}

    function GroebnerState{T1, T2, T3}(
        params
    ) where {T1 <: CoeffZZ, T2 <: CoeffQQ, T3 <: CoeffZp}
        new(
            Vector{Vector{T1}}(),
            Vector{Vector{T2}}(),
            Vector{Vector{Vector{T3}}}(),
            0,
            Vector{T1}(),
            Vector{T1}(),
            Vector{T2}(),
            Vector{BitVector}(),
            Vector{BitVector}(),
            CoefficientBuffer(),
            Vector{Vector{Vector{Vector{T3}}}}(),
            Vector{Vector{Vector{T1}}}(),
            Vector{Vector{Vector{T2}}}()
        )
    end
end

function gb_modular_select_indices0(gb_coeffs, params::AlgorithmParameters)
    indices_selection = Vector{Tuple{Int, Int}}(undef, length(gb_coeffs))
    rng = params.rng
    k = 1
    tot = 0
    @inbounds for i in 1:length(gb_coeffs)
        l = length(gb_coeffs[i])
        tot += l
        isone(l) && continue
        nl = max(floor(Int, log2(l)) - 1, 1)
        while k + nl + 1 > length(indices_selection)
            resize!(indices_selection, 2 * length(indices_selection))
        end
        for j in 1:nl
            indices_selection[k] = (i, rand(rng, 2:l))
            k += 1
        end
        indices_selection[k] = (i, l)
        k += 1
    end
    resize!(indices_selection, k - 1)
    unique!(indices_selection)
    tot, indices_selection
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

function _clear_denominators!(
    ctx::Context,
    buffer::CoefficientBuffer,
    coeffs_qq::Vector{Vector{T}}
) where {T <: CoeffQQ}
    coeffs_zz = [[BigInt(0) for _ in 1:length(c)] for c in coeffs_qq]
    den, buf = buffer.scalebuf1, buffer.scalebuf2
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

function clear_denominators!(
    ctx::Context,
    buffer::CoefficientBuffer,
    basis::Basis{T};
    deepcopy=false
) where {T <: CoeffQQ}
    coeffs_zz = _clear_denominators!(ctx, buffer, basis.coeffs)
    if deepcopy
        basis_deep_copy_with_new_coeffs(ctx, basis, coeffs_zz)
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
    coeffbuff::CoefficientBuffer,
    coeffs_zz::Vector{Vector{T1}},
    coeffs_ff::Vector{Vector{T2}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    p   = coeffbuff.reducebuf1
    buf = coeffbuff.reducebuf2
    c   = coeffbuff.reducebuf3
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
    ctx::Context,
    coeffbuff::CoefficientBuffer,
    ring::PolyRing,
    coeffs_zz::Vector{Vector{T1}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    coeffs_ff = [Vector{CoeffModular}(undef, length(c)) for c in coeffs_zz]
    ring_ff, coeffs_ff = reduce_modulo_p!(ring, coeffbuff, coeffs_zz, coeffs_ff, prime)
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    ctx::Context,
    buffer::CoefficientBuffer,
    ring::PolyRing,
    basis::Basis{T1},
    prime::T2;
    deepcopy=true
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    ring_ff, coeffs_ff = reduce_modulo_p!(ctx, buffer, ring, basis.coeffs, prime)
    new_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(ctx, basis, coeffs_ff)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_ff)
    end
    ring_ff, new_basis
end

function reduce_modulo_p_in_batch!(
    ctx::Context,
    coeffbuff::CoefficientBuffer,
    ring::PolyRing,
    basis::Basis{C},
    prime_xn::NTuple{N, T}
) where {C, N, T}
    coeffs_zz = basis.coeffs
    coeffs_ff_xn = [Vector{CompositeNumber{N, T}}(undef, length(c)) for c in coeffs_zz]

    p = coeffbuff.reducebuf1
    buf = coeffbuff.reducebuf2
    xn = map(_ -> BigInt(0), 1:N)
    c = coeffbuff.reducebuf3
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
    basis_ff_4x = basis_deep_copy_with_new_coeffs(ctx, basis, coeffs_ff_xn)

    ring_ff_4x, basis_ff_4x
end

# Resizes the state so that it has enough space to store the GB coefficients
# `gb_coeffs`
function resize_state_if_needed!(
    state::GroebnerState,
    gb_coeffs::Vector{Vector{T}}
) where {T}
    resize!(state.gb_coeffs_zz, length(gb_coeffs))
    resize!(state.gb_coeffs_qq, length(gb_coeffs))
    resize!(state.is_crt_reconstructed_mask, length(gb_coeffs))
    resize!(state.is_rational_reconstructed_mask, length(gb_coeffs))

    @inbounds for i in 1:length(gb_coeffs)
        state.gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.gb_coeffs_qq[i] = [Rational{BigInt}(1) for _ in 1:length(gb_coeffs[i])]
        state.is_crt_reconstructed_mask[i] = falses(length(gb_coeffs[i]))
        state.is_rational_reconstructed_mask[i] = falses(length(gb_coeffs[i]))
    end

    nothing
end

###
# CRT

function full_simultaneous_crt_reconstruct!(state::GroebnerState, lucky::LuckyPrimes)
    if isempty(state.gb_coeffs_zz)
        @log :misc "Using full trivial CRT reconstruction"
        coeffs_ff = state.gb_coeffs_ff_all[1]
        resize_state_if_needed!(state, coeffs_ff)
        gb_coeffs_zz = state.gb_coeffs_zz
        @inbounds for i in 1:length(coeffs_ff)
            for j in 1:length(coeffs_ff[i])
                Base.GMP.MPZ.set_ui!(gb_coeffs_zz[i][j], coeffs_ff[i][j])
            end
        end
        Base.GMP.MPZ.mul_ui!(lucky.modulo, lucky.primes[1])
        return nothing
    end
    # @invariant length(coeffs_ff) == length(state.gb_coeffs_zz)

    @log :misc "Using full CRT reconstruction"
    gb_coeffs_zz = state.gb_coeffs_zz
    is_crt_reconstructed_mask = state.is_crt_reconstructed_mask

    # Takes the lock..
    @invariant length(state.gb_coeffs_ff_all) == length(lucky.primes)

    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4
    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7
    M0 = buffer.reconstructbuf8
    MM0 = buffer.reconstructbuf9

    @inbounds for i in 1:length(gb_coeffs_zz)
        for j in 1:length(gb_coeffs_zz[i])
            if is_crt_reconstructed_mask[i][j]
                continue
            end
            Base.GMP.MPZ.set_ui!(gb_coeffs_zz[i][j], state.gb_coeffs_ff_all[1][i][j])
        end
    end

    @inbounds for i in 1:length(gb_coeffs_zz)
        if is_crt_reconstructed_mask[i][1]
            continue
        end
        Base.GMP.MPZ.set_ui!(gb_coeffs_zz[i][1], CoeffModular(1))
    end

    n = length(lucky.primes)
    rems = Vector{UInt64}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.primes
    crt_precompute!(M, n1, n2, mults, moduli)

    @log :debug "Using simultaneous CRT with moduli $moduli"

    @inbounds for i in 1:length(gb_coeffs_zz)
        for j in 2:length(gb_coeffs_zz[i])
            if is_crt_reconstructed_mask[i][j]
                continue
            end

            for ell in 1:length(lucky.primes)
                rems[ell] = state.gb_coeffs_ff_all[ell][i][j] % UInt64
            end

            crt!(M, buf, n1, n2, rems, mults)

            # cf_zz = selected_prev_coeffs_zz[i]
            # crt!(MM0, buf, n1, n2, cf_zz, invm1, buf, invm2)

            # Base.GMP.MPZ.set!(selected_coeffs_zz[i], buf)

            Base.GMP.MPZ.set!(gb_coeffs_zz[i][j], buf)
        end
    end

    Base.GMP.MPZ.set!(lucky.modulo, M)

    nothing
end

function full_simultaneous_crt_reconstruct_changematrix!(
    state::GroebnerState,
    lucky::LuckyPrimes
)
    if isempty(state.changematrix_coeffs_zz)
        @log :misc "Using full trivial CRT reconstruction"
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
        # Base.GMP.MPZ.mul_ui!(lucky.modulo, lucky.primes[1])
        return nothing
    end
    # @invariant length(coeffs_ff) == length(state.changematrix_coeffs_zz)

    @log :misc "Using full CRT reconstruction"
    changematrix_coeffs_zz = state.changematrix_coeffs_zz

    # Takes the lock..
    @invariant length(state.changematrix_coeffs_ff_all) == length(lucky.primes)

    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4
    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7
    M0 = buffer.reconstructbuf8
    MM0 = buffer.reconstructbuf9

    n = length(lucky.primes)
    rems = Vector{UInt64}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.primes
    crt_precompute!(M, n1, n2, mults, moduli)

    @log :debug "Using simultaneous CRT with moduli $moduli"

    @inbounds for i in 1:length(changematrix_coeffs_zz)
        for j in 1:length(changematrix_coeffs_zz[i])
            for k in 1:length(changematrix_coeffs_zz[i][j])
                for ell in 1:length(lucky.primes)
                    rems[ell] = state.changematrix_coeffs_ff_all[ell][i][j][k] % UInt64
                end
                crt!(M, buf, n1, n2, rems, mults)
                Base.GMP.MPZ.set!(changematrix_coeffs_zz[i][j][k], buf)
            end
        end
    end

    # Base.GMP.MPZ.set!(lucky.modulo, M)

    nothing
end

@timeit function partial_incremental_crt_reconstruct!(
    state::GroebnerState,
    lucky::LuckyPrimes,
    indices_selection::Vector{Tuple{Int, Int}}
)
    gb_coeffs_ff = state.gb_coeffs_ff_all[end]
    selected_coeffs_zz = state.selected_coeffs_zz
    selected_prev_coeffs_zz = state.selected_prev_coeffs_zz
    selected_coeffs_qq = state.selected_coeffs_qq
    gb_coeffs_zz = state.gb_coeffs_zz
    is_crt_reconstructed_mask = state.is_crt_reconstructed_mask

    if length(selected_prev_coeffs_zz) < length(indices_selection)
        @log :misc "Using partial incremental trivial CRT"
        resize!(selected_prev_coeffs_zz, length(indices_selection))
        resize!(selected_coeffs_zz, length(indices_selection))
        resize!(selected_coeffs_qq, length(indices_selection))
        @inbounds for i in 1:length(indices_selection)
            i1, i2 = indices_selection[i]
            selected_prev_coeffs_zz[i] = BigInt(0)
            selected_coeffs_zz[i] = BigInt(0)
            selected_coeffs_qq[i] = Rational{BigInt}(0)

            Base.GMP.MPZ.set_ui!(selected_coeffs_zz[i], gb_coeffs_ff[i1][i2])
            Base.GMP.MPZ.set!(gb_coeffs_zz[i1][i2], selected_coeffs_zz[i])
        end
        state.prev_index = 1
        return nothing
    end

    @log :debug "Using partial incremental CRT"

    @inbounds for i in 1:length(indices_selection)
        Base.GMP.MPZ.set!(selected_prev_coeffs_zz[i], selected_coeffs_zz[i])
    end

    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4
    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7

    crt_precompute!(M, n1, n2, invm1, lucky.modulo, invm2, last(lucky.primes))

    @inbounds for i in 1:length(indices_selection)
        i1, i2 = indices_selection[i]

        c_zz = selected_coeffs_zz[i]
        c_ff = gb_coeffs_ff[i1][i2] % UInt64
        crt!(M, buf, n1, n2, c_zz, invm1, c_ff, invm2)

        Base.GMP.MPZ.set!(selected_coeffs_zz[i], buf)

        # Mark that the coefficient is already reconstructed and save it
        is_crt_reconstructed_mask[i1][i2] = true
        Base.GMP.MPZ.set!(gb_coeffs_zz[i1][i2], selected_coeffs_zz[i])
    end

    Base.GMP.MPZ.mul_ui!(lucky.modulo, last(lucky.primes))
    state.prev_index += 1
    @log :debug "After:" lucky.modulo state.prev_index

    nothing
end

function partial_simultaneous_crt_reconstruct!(
    state::GroebnerState{T1, T2, T3},
    lucky::LuckyPrimes,
    indices_selection::Vector{Tuple{Int, Int}}
) where {T1, T2, T3}
    selected_coeffs_zz = state.selected_coeffs_zz
    selected_prev_coeffs_zz = state.selected_prev_coeffs_zz
    selected_coeffs_qq = state.selected_coeffs_qq
    is_crt_reconstructed_mask = state.is_crt_reconstructed_mask

    if length(selected_prev_coeffs_zz) < length(indices_selection)
        partial_incremental_crt_reconstruct!(state, lucky, indices_selection)
        return nothing
    end

    prev_index = state.prev_index
    n = length(lucky.primes) - prev_index
    @invariant n > 0
    if n == 1
        @log :misc "Since there is only 1 new modulo, using incremental CRT"
        partial_incremental_crt_reconstruct!(state, lucky, indices_selection)
        return nothing
    end

    @log :misc "Using partial simultaneous CRT on range $(prev_index + 1)..$(length(lucky.primes))"

    @invariant n > 1
    @inbounds for i in 1:length(indices_selection)
        Base.GMP.MPZ.set!(selected_prev_coeffs_zz[i], selected_coeffs_zz[i])
    end

    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4

    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7
    M0 = buffer.reconstructbuf8
    MM0 = buffer.reconstructbuf9

    rems = Vector{UInt}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.primes[(prev_index + 1):end]
    crt_precompute!(M, n1, n2, mults, moduli)

    Base.GMP.MPZ.set!(M0, lucky.modulo)
    crt_precompute!(MM0, n1, n2, invm1, M0, invm2, M)

    @log :debug "Using simultaneous CRT with moduli $moduli" lucky.modulo M0 M MM0

    @inbounds for i in 1:length(indices_selection)
        i1, i2 = indices_selection[i]

        for j in (prev_index + 1):length(lucky.primes)
            rems[j - prev_index] = state.gb_coeffs_ff_all[j][i1][i2] % UInt64
        end

        crt!(M, buf, n1, n2, rems, mults)

        cf_zz = selected_prev_coeffs_zz[i]
        crt!(MM0, buf, n1, n2, cf_zz, invm1, buf, invm2)

        Base.GMP.MPZ.set!(selected_coeffs_zz[i], buf)

        is_crt_reconstructed_mask[i1][i2] = true
        Base.GMP.MPZ.set!(state.gb_coeffs_zz[i1][i2], selected_coeffs_zz[i])
    end

    state.prev_index += n
    Base.GMP.MPZ.set!(lucky.modulo, MM0)

    nothing
end

# Reconstruct coefficients of the basis using rational reconstrction.
#
# state.gb_coeffs_zz -- coefficients of the basis modulo P
# state.gb_coeffs_qq -- coefficients of the basis in the rationals
# 
# Writes the coefficients of the basis modulo P reconstructed to rational
# numbers to state.gb_coeffs_qq inplace
#
# Returns true is the reconstrction is successfull, false otherwise.
@timeit function full_rational_reconstruct!(
    state::GroebnerState,
    lucky::LuckyPrimes,
    use_flint::Bool
)
    modulo = lucky.modulo
    @invariant modulo == prod(BigInt, lucky.primes)

    buffer = state.buffer
    bnd = ratrec_reconstruction_bound(modulo)

    buf, buf1 = buffer.reconstructbuf1, buffer.reconstructbuf2
    buf2, buf3 = buffer.reconstructbuf3, buffer.reconstructbuf4
    u1, u2 = buffer.reconstructbuf5, buffer.reconstructbuf6
    u3, v1 = buffer.reconstructbuf7, buffer.reconstructbuf8
    v2, v3 = buffer.reconstructbuf9, buffer.reconstructbuf10
    gb_coeffs_zz = state.gb_coeffs_zz
    gb_coeffs_qq = state.gb_coeffs_qq
    is_rational_reconstructed_mask = state.is_rational_reconstructed_mask

    @invariant length(gb_coeffs_zz) == length(gb_coeffs_qq)

    if use_flint
        nemo_bnd = Nemo.ZZRingElem(bnd)
        nemo_modulo = Nemo.ZZRingElem(modulo)

        @inbounds for i in 1:length(gb_coeffs_zz)
            @invariant length(gb_coeffs_zz[i]) == length(gb_coeffs_qq[i])
            # Skip reconstrction of the first coefficient, it is equal to one in the
            # reduced basis
            for j in 2:length(gb_coeffs_zz[i])
                if is_rational_reconstructed_mask[i][j]
                    continue
                end

                cz = gb_coeffs_zz[i][j]
                nemo_rem = Nemo.ZZRingElem(cz)
                success, pq = ratrec_nemo_2(nemo_rem, nemo_modulo, nemo_bnd, nemo_bnd)
                !success && return false

                gb_coeffs_qq[i][j] = pq
            end

            @invariant isone(gb_coeffs_qq[i][1])
        end
    else
        @inbounds for i in 1:length(gb_coeffs_zz)
            @invariant length(gb_coeffs_zz[i]) == length(gb_coeffs_qq[i])
            # Skip reconstrction of the first coefficient, it is equal to one in the
            # reduced basis
            for j in 2:length(gb_coeffs_zz[i])
                if is_rational_reconstructed_mask[i][j]
                    continue
                end

                cz = gb_coeffs_zz[i][j]
                cq = gb_coeffs_qq[i][j]
                num, den = numerator(cq), denominator(cq)
                success = ratrec!(
                    num,
                    den,
                    bnd,
                    buf,
                    buf1,
                    buf2,
                    buf3,
                    u1,
                    u2,
                    u3,
                    v1,
                    v2,
                    v3,
                    cz,
                    modulo
                )

                !success && return false
            end

            @invariant isone(gb_coeffs_qq[i][1])
        end
    end

    true
end

function full_rational_reconstruct_changematrix!(
    state::GroebnerState,
    lucky::LuckyPrimes,
    use_flint::Bool
)
    modulo = lucky.modulo
    @invariant modulo == prod(BigInt, lucky.primes)

    buffer = state.buffer
    bnd = ratrec_reconstruction_bound(modulo)

    buf, buf1 = buffer.reconstructbuf1, buffer.reconstructbuf2
    buf2, buf3 = buffer.reconstructbuf3, buffer.reconstructbuf4
    u1, u2 = buffer.reconstructbuf5, buffer.reconstructbuf6
    u3, v1 = buffer.reconstructbuf7, buffer.reconstructbuf8
    v2, v3 = buffer.reconstructbuf9, buffer.reconstructbuf10
    changematrix_coeffs_zz = state.changematrix_coeffs_zz
    changematrix_coeffs_qq = state.changematrix_coeffs_qq
    # is_rational_reconstructed_mask = state.is_rational_reconstructed_mask

    @invariant length(changematrix_coeffs_zz) == length(changematrix_coeffs_qq)
    @assert use_flint
    nemo_modulo = Nemo.ZZRingElem(modulo)
    nemo_bnd = Nemo.ZZRingElem(bnd)

    @inbounds for i in 1:length(changematrix_coeffs_zz)
        @invariant length(changematrix_coeffs_zz[i]) == length(changematrix_coeffs_qq[i])
        # !!!
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

@timeit function partial_rational_reconstruct!(
    state::GroebnerState,
    lucky::LuckyPrimes,
    indices_selection::Vector{Tuple{Int, Int}},
    use_flint::Bool
)
    modulo = lucky.modulo
    @invariant modulo == prod(BigInt, lucky.primes)

    buffer = state.buffer
    bnd = ratrec_reconstruction_bound(modulo)

    buf, buf1 = buffer.reconstructbuf1, buffer.reconstructbuf2
    buf2, buf3 = buffer.reconstructbuf3, buffer.reconstructbuf4
    u1, u2 = buffer.reconstructbuf5, buffer.reconstructbuf6
    u3, v1 = buffer.reconstructbuf7, buffer.reconstructbuf8
    v2, v3 = buffer.reconstructbuf9, buffer.reconstructbuf10

    selected_coeffs_zz = state.selected_coeffs_zz
    selected_coeffs_qq = state.selected_coeffs_qq
    gb_coeffs_qq = state.gb_coeffs_qq
    is_rational_reconstructed_mask = state.is_rational_reconstructed_mask

    if use_flint
        nemo_modulo = Nemo.ZZRingElem(modulo)
        nemo_bnd = Nemo.ZZRingElem(bnd)

        @inbounds for i in 1:length(indices_selection)
            i1, i2 = indices_selection[i]
            cz = selected_coeffs_zz[i]
            nemo_rem = Nemo.ZZRingElem(cz)

            success, pq = ratrec_nemo_2(nemo_rem, nemo_modulo, nemo_bnd, nemo_bnd)
            !success && return false

            selected_coeffs_qq[i] = pq

            # Mark that the coefficient is already reconstructed
            is_rational_reconstructed_mask[i1][i2] = true
            tnum, tden = numerator(gb_coeffs_qq[i1][i2]), denominator(gb_coeffs_qq[i1][i2])
            Base.GMP.MPZ.set!(tnum, numerator(pq))
            Base.GMP.MPZ.set!(tden, denominator(pq))
        end
    else
        @inbounds for i in 1:length(indices_selection)
            i1, i2 = indices_selection[i]

            cz = selected_coeffs_zz[i]
            cq = selected_coeffs_qq[i]
            num, den = numerator(cq), denominator(cq)
            success = ratrec!(
                num,
                den,
                bnd,
                buf,
                buf1,
                buf2,
                buf3,
                u1,
                u2,
                u3,
                v1,
                v2,
                v3,
                cz,
                modulo
            )

            !success && return false

            # Mark that the coefficient is already reconstructed
            is_rational_reconstructed_mask[i1][i2] = true
            tnum, tden = numerator(gb_coeffs_qq[i1][i2]), denominator(gb_coeffs_qq[i1][i2])
            Base.GMP.MPZ.set!(tnum, num)
            Base.GMP.MPZ.set!(tden, den)
        end
    end

    true
end
