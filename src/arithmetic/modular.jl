
#=
    The file contains functions to operate modular reconstrction and CRT
=#

#------------------------------------------------------------------------------

# Rational number reconstruction implementation borrowed from CLUE
# and modified a bit to suit the 'Modern Computer Algebra' definitions
# Returns a rational r // h of QQ field in a canonical form such that
#   r // h ≡ a (mod m)
#
# let n = max( λ(a), λ(m) ) , where λ(x) is a number of bits for x
# O(n^2)
function rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}
    a = mod(a, m)
    if a == 0 || m == 0
        return QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return QQ(1, 1)
    end
    bnd = sqrt(float(m) / 2)

    # an issue: we can not make good use of staticarrays here
    # as long as BigInts are heap allocated anyways
    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        # TODO: use MutableArithmetics
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])
    # changed from `<= bnd` to `<= m / bnd`
    # we can speed up this !
    if t <= m / bnd && gcd(r, t) == 1
        return QQ(r, t)
    end

    # TODO: not needed
    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))

    return QQ(0, 1)
end

function reconstruct_modulo(
        coeffs_zz::Vector{Vector{BigInt}},
        modulo::BigInt)

    coeffs_qq = Vector{Vector{Rational{BigInt}}}(undef, length(coeffs_zz))
    for i in 1:length(coeffs_zz)
        poly_zz = coeffs_zz[i]
        coeffs_qq[i] = Vector{Rational{BigInt}}(undef, length(poly_zz))
        for j in 1:length(poly_zz)
            coeffs_qq[i][j] = rational_reconstruction(poly_zz[j], modulo)
        end
    end

    coeffs_qq
end

#------------------------------------------------------------------------------

function CRT(
            a1::BigInt, m1::BigInt,
            a2::U, m2::BigInt) where {U<:Integer}

    CRT(a1, m1, BigInt(a2), m2)
end

function CRT(a1::BigInt, m1::BigInt, a2::BigInt, m2::BigInt)
    M = m1 * m2
    Ms1, Ms2 = m2, m1
    t1 = invmod(Ms1, m1)
    t2 = invmod(Ms2, m2)
    mod(a1*t1*Ms1 + a2*t2*Ms2, M)::BigInt
end

# TODO: move to coeffs.jl
function reconstruct_crt!(
            gbcoeffs_ff::Vector{Vector{UInt64}}, ch::UInt64)

    gbcoeffs_zz = Vector{Vector{BigInt}}(undef, length(gbcoeffs_ff))
    bigch = BigInt(ch)
    for i in 1:length(gbcoeffs_ff)
        poly_ff = gbcoeffs_ff[i]
        gbcoeffs_zz[i] = Vector{BigInt}(undef, length(poly_ff))
        for j in 1:length(poly_ff)
            gbcoeffs_zz[i][j] = BigInt(poly_ff[j])
        end
    end
    gbcoeffs_zz, bigch
end

function reconstruct_crt!(
            gbcoeffs_accum::Vector{Vector{BigInt}}, modulo::BigInt,
            gbcoeffs_ff::Vector{Vector{UInt64}}, ch::UInt64)

    if isempty(gbcoeffs_accum)
        return reconstruct_crt!(gbcoeffs_ff, ch)
    end

    gbcoeffs_zz = Vector{Vector{BigInt}}(undef, length(gbcoeffs_ff))
    bigch = BigInt(ch)

    for i in 1:length(gbcoeffs_ff)
        poly_zz = gbcoeffs_accum[i]
        poly_ff = gbcoeffs_ff[i]
        gbcoeffs_zz[i] = Vector{BigInt}(undef, length(poly_ff))
        for j in 1:length(poly_ff)
            gbcoeffs_zz[i][j] = CRT(poly_zz[j], modulo, poly_ff[j], bigch)
        end
    end
    # TODO: MutableArithmetics
    gbcoeffs_zz, modulo * bigch
end

#------------------------------------------------------------------------------
