
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
mutable struct GroebnerState{T1 <: CoeffZZ, T2 <: CoeffQQ}
    gb_coeffs_zz::Vector{Vector{T1}}
    prev_gb_coeffs_zz::Vector{Vector{T1}}
    gb_coeffs_qq::Vector{Vector{T2}}
    buffer::CoefficientBuffer
    function GroebnerState{T1, T2}(params) where {T1 <: CoeffZZ, T2 <: CoeffQQ}
        new(
            Vector{Vector{T1}}(undef, 0),
            Vector{Vector{T1}}(undef, 0),
            Vector{Vector{T2}}(undef, 0),
            CoefficientBuffer()
        )
    end
end

function majority_vote!(state, basis_ff, tracer, params)
    true
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

# TODO: scale numerators inplace!
function clear_denominators!(
    buffer::CoefficientBuffer,
    coeffs_zz::Vector{Vector{C1}},
    coeffs_qq::Vector{Vector{C2}}
) where {C1 <: CoeffZZ, C2 <: CoeffQQ}
    @assert length(coeffs_zz) == length(coeffs_qq)
    den, buf = buffer.scalebuf1, buffer.scalebuf2
    @inbounds for i in 1:length(coeffs_qq)
        @assert length(coeffs_zz[i]) == length(coeffs_qq[i])
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
    copy_basis(basis, coeffs_zz, deepcopy=deepcopy)
end

function reduce_modulo_p!(
    ring,
    coeffbuff::CoefficientBuffer,
    coeffs_zz::Vector{Vector{T1}},
    coeffs_ff::Vector{Vector{T2}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffFF}
    p   = coeffbuff.reducebuf1
    buf = coeffbuff.reducebuf2
    c   = coeffbuff.reducebuf3
    Base.GMP.MPZ.set_ui!(p, prime)

    for i in 1:length(coeffs_zz)
        cfs_zz_i = coeffs_zz[i]
        for j in 1:length(cfs_zz_i)
            Base.GMP.MPZ.set!(c, cfs_zz_i[j])
            if Base.GMP.MPZ.cmp_ui(c, 0) < 0
                Base.GMP.MPZ.fdiv_q!(buf, c, p)
                Base.GMP.MPZ.neg!(buf)
                Base.GMP.MPZ.mul_ui!(buf, prime)
                Base.GMP.MPZ.add!(c, buf)
            end
            Base.GMP.MPZ.tdiv_r!(buf, c, p)
            coeffs_ff[i][j] = CoeffModular(buf)
        end
    end
    ring_ff = PolyRing(ring.nvars, ring.ord, Int(prime))
    ring_ff, coeffs_ff
end

function reduce_modulo_p!(
    coeffbuff::CoefficientBuffer,
    ring,
    coeffs_zz::Vector{Vector{T1}},
    prime::T2
) where {T1 <: CoeffZZ, T2 <: CoeffFF}
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
    ring_ff, copy_basis(basis, coeffs_ff, deepcopy=deepcopy)
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
    @inbounds for i in 1:length(gb_coeffs)
        # TODO: we are being too generous with memory here
        state.gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.prev_gb_coeffs_zz[i] = [BigInt(0) for _ in 1:length(gb_coeffs[i])]
        state.gb_coeffs_qq[i] = [Rational{BigInt}(1) for _ in 1:length(gb_coeffs[i])]
    end
    nothing
end

# Reconstruct coefficients of the basis using CRT.
#
# state.gb_coeffs_zz -- coefficients of the basis modulo P1*P2*...*Pk.
# basis_ff.coeffs -- coefficients of the basis modulo new prime P.
# 
# Writes the coefficients of the basis modulo P * P1*P2*...*Pk to
# state.gb_coeffs_zz inplace
function crt_reconstruct!(
    state::GroebnerState,
    ring::PolyRing,
    lucky::LuckyPrimes,
    basis_ff::Basis{T}
) where {T <: CoeffFF}
    if isempty(state.gb_coeffs_zz)
        # If first time reconstruction
        resize_state_if_needed!(state, basis_ff.coeffs)
        crt_reconstruct_trivial!(state, basis_ff.coeffs)
        Base.GMP.MPZ.mul_ui!(lucky.modulo, last(lucky.primes))
        return nothing
    end
    @invariant length(basis_ff.coeffs) == length(state.gb_coeffs_zz)
    gb_coeffs_zz = state.gb_coeffs_zz
    prev_gb_coeffs_zz = state.prev_gb_coeffs_zz
    # Copy current basis coefficients to the previous array
    @inbounds for i in 1:length(gb_coeffs_zz)
        @invariant length(gb_coeffs_zz[i]) == length(prev_gb_coeffs_zz[i])
        for j in 1:length(gb_coeffs_zz[i])
            Base.GMP.MPZ.set!(prev_gb_coeffs_zz[i][j], gb_coeffs_zz[i][j])
        end
    end
    # Set the buffers for CRT and precompute some values
    buffer = state.buffer
    buf = buffer.reconstructbuf1
    n1, n2 = buffer.reconstructbuf2, buffer.reconstructbuf3
    M = buffer.reconstructbuf4
    bigch = buffer.reconstructbuf5
    invm1, invm2 = buffer.reconstructbuf6, buffer.reconstructbuf7
    characteristic = ring.ch
    Base.GMP.MPZ.set_ui!(bigch, characteristic)
    Base.GMP.MPZ.mul_ui!(M, lucky.modulo, characteristic)
    Base.GMP.MPZ.gcdext!(buf, invm1, invm2, lucky.modulo, bigch)
    gb_coeffs_ff = basis_ff.coeffs
    @invariant length(gb_coeffs_zz) == length(gb_coeffs_ff)
    @inbounds for i in 1:length(gb_coeffs_zz)
        @invariant length(gb_coeffs_zz[i]) == length(gb_coeffs_ff[i])
        for j in 1:length(gb_coeffs_zz[i])
            ca = gb_coeffs_zz[i][j]
            cf = gb_coeffs_ff[i][j]
            CRT!(M, buf, n1, n2, ca, invm1, cf, invm2, lucky.modulo, bigch)
            Base.GMP.MPZ.set!(gb_coeffs_zz[i][j], buf)
        end
    end
    Base.GMP.MPZ.mul_ui!(lucky.modulo, last(lucky.primes))
    nothing
end

function crt_reconstruct_trivial!(
    state::GroebnerState,
    gb_coeffs_ff::Vector{Vector{T1}}
) where {T1 <: CoeffFF}
    gb_coeffs_zz = state.gb_coeffs_zz
    @inbounds for i in 1:length(gb_coeffs_ff)
        for j in 1:length(gb_coeffs_ff[i])
            Base.GMP.MPZ.set_ui!(gb_coeffs_zz[i][j], gb_coeffs_ff[i][j])
        end
    end
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
function rational_reconstruct!(state::GroebnerState, lucky::LuckyPrimes)
    modulo = lucky.modulo
    buffer     = state.buffer
    bnd        = rational_reconstruction_bound(modulo)
    buf, buf1  = buffer.reconstructbuf1, buffer.reconstructbuf2
    buf2, buf3 = buffer.reconstructbuf3, buffer.reconstructbuf4
    u1, u2     = buffer.reconstructbuf5, buffer.reconstructbuf6
    u3, v1     = buffer.reconstructbuf7, buffer.reconstructbuf8
    v2, v3     = buffer.reconstructbuf9, buffer.reconstructbuf10
    gb_coeffs_zz = state.gb_coeffs_zz
    gb_coeffs_qq = state.gb_coeffs_qq
    @invariant length(gb_coeffs_zz) == length(gb_coeffs_qq)
    @inbounds for i in 1:length(gb_coeffs_zz)
        @invariant length(gb_coeffs_zz[i]) == length(gb_coeffs_qq[i])
        # skip reconstrction of the first coefficient, it is equal to one in the
        # reduced basis
        for j in 2:length(gb_coeffs_zz[i])
            cz = gb_coeffs_zz[i][j]
            cq = gb_coeffs_qq[i][j]
            num, den = numerator(cq), denominator(cq)
            success = rational_reconstruction!(
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
    true
end
