
#=
    `CoeffBuffer` contains BigInt buffers used to minimize memory allocations
    during coefficient scaling, rational reconstruction and CRT reconstruction
=#
mutable struct CoeffBuffer
    # wow, 40 bytes per one BigInt

    # buffers that should be used
    # only for coefficient scaling
    scalebuf1::BigInt
    scalebuf2::BigInt

    # buffers that should be used
    # only for coefficient reduction
    reducebuf1::BigInt
    reducebuf2::BigInt
    reducebuf3::BigInt

    # buffers that should be used
    # only for coefficient reconstrction
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

    function CoeffBuffer()
        new(BigInt(), BigInt(), BigInt(), BigInt(), BigInt(),
            BigInt(), BigInt(), BigInt(), BigInt(), BigInt(), BigInt(),
            BigInt(), BigInt(), BigInt(), BigInt())
    end
end

#------------------------------------------------------------------------------

#=
    The structure stores accumulated integer coefficients and
    accumulated rational coefficients of a groebner basis
=#
mutable struct CoeffAccum
    gb_coeffs_zz::Vector{Vector{CoeffZZ}}
    prev_gb_coeffs_zz::Vector{Vector{CoeffZZ}}
    gb_coeffs_qq::Vector{Vector{CoeffQQ}}

    function CoeffAccum()
        new(Vector{Vector{CoeffZZ}}(undef, 0),
            Vector{Vector{CoeffZZ}}(undef, 0),
            Vector{Vector{CoeffQQ}}(undef, 0))
    end
end

#------------------------------------------------------------------------------

# compute gcd of denominators of `coeffs`
function common_denominator(coeffbuff::CoeffBuffer, coeffs::Vector{CoeffQQ})
    den = coeffbuff.scalebuf1
    Base.GMP.MPZ.set_si!(den, 1)
    for c in coeffs
        # TODO
        # use buffer?
        Base.GMP.MPZ.lcm!(den, denominator(c))
    end
    den
end

# compute gcd of denominators of `coeffs`
function common_denominator(coeffs::Vector{CoeffQQ})
    den = BigInt()
    Base.GMP.MPZ.set_si!(den, 1)
    for c in coeffs
        # TODO
        # use buffer?
        Base.GMP.MPZ.lcm!(den, denominator(c))
    end
    den
end

# scales denominators inplace
function scale_denominators!(
            coeffbuff::CoeffBuffer,
            coeffs_qq::Vector{Vector{CoeffQQ}},
            coeffs_zz::Vector{Vector{CoeffZZ}})

    @assert length(coeffs_zz) == length(coeffs_qq)

    buf = coeffbuff.scalebuf2
    for i in 1:length(coeffs_qq)
        @assert length(coeffs_zz[i]) == length(coeffs_qq[i])

        den = common_denominator(coeffbuff, coeffs_qq[i])
        sz  = Base.GMP.MPZ.sizeinbase(den, 2)
        @inbounds for j in 1:length(coeffs_qq[i])
            num = numerator(coeffs_qq[i][j])
            Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs_qq[i][j]))
            Base.GMP.MPZ.realloc2!(coeffs_zz[i][j], sz)
            Base.GMP.MPZ.mul!(coeffs_zz[i][j], num, buf)
        end
    end

    coeffs_zz
end

function scale_denominators(
            coeffbuff::CoeffBuffer,
            coeffs_qq::Vector{Vector{CoeffQQ}})
    coeffs_zz = [[CoeffZZ(0) for _ in 1:length(c)] for c in coeffs_qq]
    scale_denominators!(coeffbuff, coeffs_qq, coeffs_zz)
end

function scale_denominators(coeffs_qq::Vector{Vector{CoeffQQ}})
    coeffs_zz = [[CoeffZZ(0) for _ in 1:length(c)] for c in coeffs_qq]
    buf = BigInt()
    for i in 1:length(coeffs_qq)
        @assert length(coeffs_zz[i]) == length(coeffs_qq[i])

        den = common_denominator(coeffs_qq[i])
        @inbounds for j in 1:length(coeffs_qq[i])
            num = numerator(coeffs_qq[i][j])
            Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs_qq[i][j]))
            Base.GMP.MPZ.mul!(coeffs_zz[i][j], num, buf)
        end
    end

    coeffs_zz
end

function scale_denominators(coeffs_qq::Vector{CoeffQQ})
    coeffs_zz = [CoeffZZ(0) for _ in 1:length(coeffs_qq)]
    buf = BigInt()
    den = common_denominator(coeffs_qq)
    for i in 1:length(coeffs_qq)
        num = numerator(coeffs_qq[i])
        Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs_qq[i]))
        Base.GMP.MPZ.mul!(coeffs_zz[i], num, buf)
    end
    coeffs_zz
end

#------------------------------------------------------------------------------

function reduce_modulo!(
        coeffbuff::CoeffBuffer,
        coeffs_zz::Vector{Vector{CoeffZZ}},
        coeffs_ff::Vector{Vector{CoeffFF}},
        prime::CoeffFF)

    p   = coeffbuff.reducebuf1
    buf = coeffbuff.reducebuf2
    c   = coeffbuff.reducebuf3
    Base.GMP.MPZ.set_ui!(p, prime)

    for i in 1:length(coeffs_zz)
        cfs_zz_i = coeffs_zz[i]
        @inbounds for j in 1:length(cfs_zz_i)
            Base.GMP.MPZ.set!(c, cfs_zz_i[j])
            # TODO: rewrite this
            if Base.GMP.MPZ.cmp_ui(c, 0) < 0
                Base.GMP.MPZ.fdiv_q!(buf, c, p)
                Base.GMP.MPZ.neg!(buf)
                Base.GMP.MPZ.mul_ui!(buf, prime)
                Base.GMP.MPZ.add!(c, buf)
            end
            # @assert c >= 0 # TODO
            Base.GMP.MPZ.tdiv_r!(buf, c, p)
            coeffs_ff[i][j] = UInt64(buf)
        end
    end
end

function reduce_modulo(
        coeffbuff::CoeffBuffer,
        coeffs_zz::Vector{Vector{CoeffZZ}},
        prime::UInt64)
    coeffs_ff =  [Vector{CoeffFF}(undef, length(c)) for c in coeffs_zz]
    reduce_modulo!(coeffbuff, coeffs_zz, coeffs_ff, prime)
end

#------------------------------------------------------------------------------

function resize_accum!(coeffaccum::CoeffAccum, gb_coeffs)
    resize!(coeffaccum.gb_coeffs_zz, length(gb_coeffs))
    resize!(coeffaccum.prev_gb_coeffs_zz, length(gb_coeffs))
    resize!(coeffaccum.gb_coeffs_qq, length(gb_coeffs))
    @inbounds for i in 1:length(gb_coeffs)
        coeffaccum.gb_coeffs_zz[i] = [CoeffZZ(0) for _ in 1:length(gb_coeffs[i])]
        coeffaccum.prev_gb_coeffs_zz[i] = [CoeffZZ(0) for _ in 1:length(gb_coeffs[i])]
        coeffaccum.gb_coeffs_qq[i] = [CoeffQQ(1) for _ in 1:length(gb_coeffs[i])]
    end
end

function reconstruct_trivial_crt!(coeffbuff::CoeffBuffer,
                                    coeffaccum::CoeffAccum,
                                    gb_coeffs_ff::Vector{Vector{CoeffFF}})
    gb_coeffs_zz = coeffaccum.gb_coeffs_zz
    for i in 1:length(gb_coeffs_ff)
        @inbounds for j in 1:length(gb_coeffs_ff[i])
            Base.GMP.MPZ.set_ui!(gb_coeffs_zz[i][j], gb_coeffs_ff[i][j])
        end
    end
end

function assure_structure(coeffaccum, gb_coeffs_ff)
    if length(coeffaccum.gb_coeffs_zz) != gb_coeffs_ff
        return false
    end
    @inbounds for i in 1:length(gb_coeffs_ff)
        if length(coeffaccum.gb_coeffs_zz[i]) != gb_coeffs_ff[i]
            return false
        end
    end
    return true
end

function reconstruct_crt!(
        coeffbuff::CoeffBuffer,
        coeffaccum::CoeffAccum,
        primetracker::PrimeTracker,
        gb_coeffs_ff::Vector{Vector{CoeffFF}},
        ch::UInt64)

    if isempty(coeffaccum.gb_coeffs_zz)
        # if first time reconstruction
        resize_accum!(coeffaccum, gb_coeffs_ff)
        reconstruct_trivial_crt!(coeffbuff, coeffaccum, gb_coeffs_ff)
    else
        gb_coeffs_zz = coeffaccum.gb_coeffs_zz
        prev_gb_coeffs_zz = coeffaccum.prev_gb_coeffs_zz

        # copy to previous gb coeffs
        for i in 1:length(gb_coeffs_zz)
            for j in 1:length(gb_coeffs_zz[i])
                Base.GMP.MPZ.set!(prev_gb_coeffs_zz[i][j], gb_coeffs_zz[i][j])
            end
        end

        buf = coeffbuff.reconstructbuf1
        n1, n2 = coeffbuff.reconstructbuf2, coeffbuff.reconstructbuf3
        M = coeffbuff.reconstructbuf4
        bigch = coeffbuff.reconstructbuf5
        invm1, invm2 = coeffbuff.reconstructbuf6, coeffbuff.reconstructbuf7

        Base.GMP.MPZ.set_ui!(bigch, ch)
        Base.GMP.MPZ.mul_ui!(M, primetracker.modulo, ch)
        Base.GMP.MPZ.gcdext!(buf, invm1, invm2, primetracker.modulo, bigch)

        for i in 1:length(gb_coeffs_ff)
            @inbounds for j in 1:length(gb_coeffs_ff[i])
                ca = gb_coeffs_zz[i][j]
                cf = gb_coeffs_ff[i][j]
                CRT!(M, buf, n1, n2, ca, invm1, cf, invm2, primetracker.modulo, bigch)
                # TODO: faster set!??
                Base.GMP.MPZ.set!(gb_coeffs_zz[i][j], buf)
            end
        end
    end
    updatemodulo!(primetracker)
end

#------------------------------------------------------------------------------

function reconstruct_modulo!(
        coeffbuff::CoeffBuffer,
        coeffaccum::CoeffAccum,
        primetracker::PrimeTracker)

    modulo = primetracker.modulo

    bnd = rational_reconstruction_bound(modulo)
    buf, buf1  = coeffbuff.reconstructbuf1, coeffbuff.reconstructbuf2
    buf2, buf3 = coeffbuff.reconstructbuf3, coeffbuff.reconstructbuf4
    u1, u2     = coeffbuff.reconstructbuf5, coeffbuff.reconstructbuf6
    u3, v1     = coeffbuff.reconstructbuf7, coeffbuff.reconstructbuf8
    v2, v3     = coeffbuff.reconstructbuf9, coeffbuff.reconstructbuf10

    gb_coeffs_zz = coeffaccum.gb_coeffs_zz
    gb_coeffs_qq = coeffaccum.gb_coeffs_qq

    for i in 1:length(gb_coeffs_zz)
        @inbounds for j in 2:length(gb_coeffs_zz[i])
            cz = gb_coeffs_zz[i][j]
            cq = gb_coeffs_qq[i][j]
            num, den = numerator(cq), denominator(cq)
            success = rational_reconstruction!(num, den, bnd, buf,
                                        buf1, buf2, buf3,
                                        u1, u2, u3, v1, v2, v3,
                                        cz, modulo)
            if !success
                # @error "not success"
                return false
            end
            # TODO
            # @assert gcd(numerator(cq), denominator(cq)) == 1
        end
    end

    return true
end
