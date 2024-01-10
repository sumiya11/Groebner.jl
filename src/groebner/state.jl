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
    prev_gb_coeffs_zz::Vector{Vector{T1}}
    gb_coeffs_qq::Vector{Vector{T2}}

    #
    gb_coeffs_ff_all::Vector{Vector{Vector{T3}}}
    prev_index::Int

    #
    selected_coeffs_zz::Vector{T1}
    selected_prev_coeffs_zz::Vector{T1}
    selected_coeffs_qq::Vector{T2}
    is_crt_reconstructed_mask::Vector{BitVector}
    is_rational_reconstructed_mask::Vector{BitVector}

    buffer::CoefficientBuffer
    function GroebnerState{T1, T2, T3}(
        params
    ) where {T1 <: CoeffZZ, T2 <: CoeffQQ, T3 <: CoeffZp}
        new(
            Vector{Vector{T1}}(),
            Vector{Vector{T1}}(),
            Vector{Vector{T2}}(),
            Vector{Vector{Vector{T3}}}(),
            0,
            Vector{T1}(),
            Vector{T1}(),
            Vector{T2}(),
            Vector{BitVector}(),
            Vector{BitVector}(),
            CoefficientBuffer()
        )
    end
end

function common_denominator!(den::BigInt, coeffs::Vector{T}) where {T <: CoeffQQ}
    Base.GMP.MPZ.set_si!(den, 1)
    for c in coeffs
        Base.GMP.MPZ.lcm!(den, denominator(c))
    end
    den
end
function common_denominator(coeffs::Vector{T}) where {T <: CoeffQQ}
    common_denominator!(BigInt(), coeffs)
end

# TODO: scale numerators inplace and do not allocate new GMP instances
function clear_denominators!(
    buffer::CoefficientBuffer,
    coeffs_zz::Vector{Vector{C1}},
    coeffs_qq::Vector{Vector{C2}}
) where {C1 <: CoeffZZ, C2 <: CoeffQQ}
    @invariant length(coeffs_zz) == length(coeffs_qq)
    den, buf = buffer.scalebuf1, buffer.scalebuf2
    @inbounds for i in 1:length(coeffs_qq)
        @invariant length(coeffs_zz[i]) == length(coeffs_qq[i])
        den = common_denominator!(den, coeffs_qq[i])
        sz  = Base.GMP.MPZ.sizeinbase(den, 2)
        for j in 1:length(coeffs_qq[i])
            num = numerator(coeffs_qq[i][j])
            Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs_qq[i][j]))
            Base.GMP.MPZ.realloc2!(coeffs_zz[i][j], sz)
            Base.GMP.MPZ.mul!(coeffs_zz[i][j], num, buf)
        end
    end
    coeffs_zz
end

function clear_denominators!(
    buffer::CoefficientBuffer,
    coeffs_qq::Vector{Vector{T}}
) where {T <: CoeffQQ}
    coeffs_zz = [[BigInt(0) for _ in 1:length(c)] for c in coeffs_qq]
    clear_denominators!(buffer, coeffs_zz, coeffs_qq)
end

function clear_denominators(coeffs_qq::Vector{Vector{T}}) where {T <: CoeffQQ}
    clear_denominators!(CoefficientBuffer(), coeffs_qq)
end

function clear_denominators(coeffs_qq::Vector{T}) where {T <: CoeffQQ}
    first(clear_denominators([coeffs_qq]))
end

function clear_denominators!(
    buffer::CoefficientBuffer,
    basis::Basis{T};
    deepcopy=false
) where {T <: CoeffQQ}
    coeffs_zz = clear_denominators!(buffer, basis.coeffs)
    if deepcopy
        basis_deep_copy_with_new_coeffs(basis, coeffs_zz)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_zz)
    end
end

# Reduce x modulo a prime
function reduce_mod_p!(buf::BigInt, x::BigInt, prime::Unsigned, prime_big::BigInt)
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
            reduce_mod_p!(buf, c, prime, p)
            coeffs_ff[i][j] = CoeffModular(buf)
        end
    end
    ring_ff = PolyRing(ring.nvars, ring.ord, UInt(prime))
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    coeffbuff::CoefficientBuffer,
    ring,
    coeffs_zz::Vector{Vector{T1}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffZp}
    coeffs_ff = [Vector{CoeffModular}(undef, length(c)) for c in coeffs_zz]
    ring_ff, coeffs_ff = reduce_modulo_p!(ring, coeffbuff, coeffs_zz, coeffs_ff, prime)
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    buffer::CoefficientBuffer,
    ring::PolyRing,
    basis::Basis{<:CoeffZZ},
    prime;
    deepcopy=true
)
    ring_ff, coeffs_ff = reduce_modulo_p!(buffer, ring, basis.coeffs, prime)
    new_basis = if deepcopy
        basis_deep_copy_with_new_coeffs(basis, coeffs_ff)
    else
        basis_shallow_copy_with_new_coeffs(basis, coeffs_ff)
    end
    ring_ff, new_basis
end

function reduce_modulo_p_in_batch!(
    coeffbuff::CoefficientBuffer,
    ring::PolyRing,
    basis::Basis{C},
    prime_xn::NTuple{N, T}
) where {C, N, T}
    coeffs_zz = basis.coeffs
    coeffs_ff_xn = [Vector{CompositeInt{N, T}}(undef, length(c)) for c in coeffs_zz]

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
                k -> T(reduce_mod_p!(buf, xn[k], UInt(prime_xn[k]), prime_big_xn[k])),
                N
            )
            coeffs_ff_xn[i][j] = CompositeInt{N, T}(data)
        end
    end
    ring_ff_4x = PolyRing(ring.nvars, ring.ord, CompositeInt{N, T}(prime_xn))
    basis_ff_4x = basis_deep_copy_with_new_coeffs(basis, coeffs_ff_xn)

    ring_ff_4x, basis_ff_4x
end

# Resizes the state so that it has enough space to store the GB coefficients
# `gb_coeffs`
function resize_state_if_needed!(
    state::GroebnerState,
    gb_coeffs::Vector{Vector{T}}
) where {T}
    resize!(state.gb_coeffs_zz, length(gb_coeffs))
    resize!(state.prev_gb_coeffs_zz, length(gb_coeffs))
    resize!(state.gb_coeffs_qq, length(gb_coeffs))
    resize!(state.is_crt_reconstructed_mask, length(gb_coeffs))
    resize!(state.is_rational_reconstructed_mask, length(gb_coeffs))

    @inbounds for i in 1:length(gb_coeffs)
        state.gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.prev_gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.gb_coeffs_qq[i] = [Rational{BigInt}(1) for _ in 1:length(gb_coeffs[i])]
        state.is_crt_reconstructed_mask[i] = falses(length(gb_coeffs[i]))
        state.is_rational_reconstructed_mask[i] = falses(length(gb_coeffs[i]))
    end

    nothing
end

###
# CRT

# Reconstruct coefficients of the basis using CRT.
#
# state.gb_coeffs_zz -- coefficients of the basis modulo P1*P2*...*Pk.
# basis_ff.coeffs -- coefficients of the basis modulo new prime P.
# 
# Writes the coefficients of the basis modulo P * P1*P2*...*Pk to
# state.gb_coeffs_zz inplace
function full_incremental_crt_reconstruct!(state::GroebnerState, lucky::LuckyPrimes)
    if isempty(state.gb_coeffs_zz)
        @log level = -2 "Using full trivial CRT reconstruction"
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

    @log level = -2 "Using full CRT reconstruction"
    gb_coeffs_zz = state.gb_coeffs_zz
    prev_gb_coeffs_zz = state.prev_gb_coeffs_zz
    is_crt_reconstructed_mask = state.is_crt_reconstructed_mask

    # Takes the lock..
    @invariant length(state.gb_coeffs_ff_all) == length(lucky.primes)

    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4
    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7

    Base.GMP.MPZ.set_ui!(lucky.modulo, lucky.primes[1])

    @inbounds for i in 1:length(gb_coeffs_zz)
        @invariant length(gb_coeffs_zz[i]) == length(prev_gb_coeffs_zz[i])
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

    for idx in 2:length(lucky.primes)
        crt_precompute!(M, n1, n2, invm1, lucky.modulo, invm2, lucky.primes[idx])
        gb_coeffs_ff = state.gb_coeffs_ff_all[idx]

        @invariant length(gb_coeffs_zz) == length(gb_coeffs_ff)
        @inbounds for i in 1:length(gb_coeffs_zz)
            @invariant length(gb_coeffs_zz[i]) == length(gb_coeffs_ff[i])

            # Skip reconstruction of the first coefficient
            for j in 2:length(gb_coeffs_zz[i])
                if is_crt_reconstructed_mask[i][j]
                    continue
                end

                # Copy current basis coefficients to the previous array
                # Base.GMP.MPZ.set!(prev_gb_coeffs_zz[i][j], gb_coeffs_zz[i][j])

                c_zz = gb_coeffs_zz[i][j]
                c_ff = gb_coeffs_ff[i][j]

                crt!(M, buf, n1, n2, c_zz, invm1, c_ff, invm2)

                Base.GMP.MPZ.set!(gb_coeffs_zz[i][j], buf)
            end
        end

        Base.GMP.MPZ.mul_ui!(lucky.modulo, lucky.primes[idx])
    end

    nothing
end

function full_simultaneous_crt_reconstruct!(state::GroebnerState, lucky::LuckyPrimes)
    if isempty(state.gb_coeffs_zz)
        @log level = -2 "Using full trivial CRT reconstruction"
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

    @log level = -2 "Using full CRT reconstruction"
    gb_coeffs_zz = state.gb_coeffs_zz
    prev_gb_coeffs_zz = state.prev_gb_coeffs_zz
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
        @invariant length(gb_coeffs_zz[i]) == length(prev_gb_coeffs_zz[i])

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
    rems = Vector{CoeffModular}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.primes
    crt_precompute!(M, n1, n2, mults, moduli)

    @log level = -6 "Using simultaneous CRT with moduli $moduli"

    @inbounds for i in 1:length(gb_coeffs_zz)
        for j in 2:length(gb_coeffs_zz[i])
            if is_crt_reconstructed_mask[i][j]
                continue
            end

            for ell in 1:length(lucky.primes)
                rems[ell] = state.gb_coeffs_ff_all[ell][i][j]
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
        @log level = -2 "Using partial incremental trivial CRT"
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

    @log level = -3 "Using partial incremental CRT"

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
        c_ff = gb_coeffs_ff[i1][i2]
        crt!(M, buf, n1, n2, c_zz, invm1, c_ff, invm2)

        Base.GMP.MPZ.set!(selected_coeffs_zz[i], buf)

        # Mark that the coefficient is already reconstructed and save it
        is_crt_reconstructed_mask[i1][i2] = true
        Base.GMP.MPZ.set!(gb_coeffs_zz[i1][i2], selected_coeffs_zz[i])
    end

    Base.GMP.MPZ.mul_ui!(lucky.modulo, last(lucky.primes))
    state.prev_index += 1
    @log level = -3 "After:" lucky.modulo state.prev_index

    nothing
end

@timeit function partial_simultaneous_crt_reconstruct!(
    state::GroebnerState,
    lucky::LuckyPrimes,
    indices_selection::Vector{Tuple{Int, Int}}
)
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
    @assert n > 0
    if n == 1
        @log level = -2 "Since there is only 1 new modulo, using incremental CRT"
        partial_incremental_crt_reconstruct!(state, lucky, indices_selection)
        return nothing
    end

    @log level = -2 "Using partial simultaneous CRT on range $(prev_index + 1)..$(length(lucky.primes))"

    @assert n > 1
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

    rems = Vector{CoeffModular}(undef, n)
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    moduli = lucky.primes[(prev_index + 1):end]
    crt_precompute!(M, n1, n2, mults, moduli)

    Base.GMP.MPZ.set!(M0, lucky.modulo)
    crt_precompute!(MM0, n1, n2, invm1, M0, invm2, M)

    @log level = -6 "Using simultaneous CRT with moduli $moduli" lucky.modulo M0 M MM0

    @inbounds for i in 1:length(indices_selection)
        i1, i2 = indices_selection[i]

        for j in (prev_index + 1):length(lucky.primes)
            rems[j - prev_index] = state.gb_coeffs_ff_all[j][i1][i2]
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
                success, (num, den) = ratrec_nemo(nemo_rem, nemo_modulo)
                gb_coeffs_qq[i][j] = Base.unsafe_rational(num, den)

                !success && return false
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

        @inbounds for i in 1:length(indices_selection)
            i1, i2 = indices_selection[i]
            cz = selected_coeffs_zz[i]
            nemo_rem = Nemo.ZZRingElem(cz)

            success, (num, den) = ratrec_nemo(nemo_rem, nemo_modulo)
            selected_coeffs_qq[i] = Base.unsafe_rational(num, den)

            !success && return false

            # Mark that the coefficient is already reconstructed
            is_rational_reconstructed_mask[i1][i2] = true
            tnum, tden = numerator(gb_coeffs_qq[i1][i2]), denominator(gb_coeffs_qq[i1][i2])
            Base.GMP.MPZ.set!(tnum, num)
            Base.GMP.MPZ.set!(tden, den)
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
