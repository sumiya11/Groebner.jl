
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

    U = I[1, 0, m]
    V = I[0, 1, a]
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
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

    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))

    return QQ(0, 1)
end

#------------------------------------------------------------------------------

function reduce_modulo(fs, modulo)
    ground = GF(modulo)
    Rin = parent(first(fs))
    Rgf, _ = PolynomialRing(ground, string.(gens(Rin)), ordering=ordering(Rin))
    ans = zeros(Rgf, length(fs))
    for i in 1:length(fs)
        ans[i] = map_coefficients(c -> ground(c), fs[i])
    end
    ans
end

#------------------------------------------------------------------------------

scale_denominators(f::MPoly) = scale_denominators([f])

function scale_denominators(fs::Array)
    Rin = parent(first(fs))
    Rzz, _ = PolynomialRing(ZZ, string.(gens(Rin)), ordering=ordering(Rin))
    ans = zeros(Rzz, length(fs))
    for i in 1:length(fs)
        scaled = map_coefficients(c -> c // content(fs[i]), fs[i])
        ans[i] = change_base_ring(ZZ, scaled, parent=Rzz)
    end
    ans
end

#------------------------------------------------------------------------------

function crt(rs::T, ms::T) where {T<:AbstractArray}
    n = length(rs)
    r1 = Nemo.ZZ(rs[1])
    m1 = Nemo.ZZ(ms[1])

    crt(r1, m1, view(rs, 2:n), view(ms, 2:n))
end

function AbstractAlgebra.crt(r1, m1, rs::T, ms::T) where {T<:AbstractArray}
    n = length(rs)
    if n == 0
        return r1
    end
    return crt(
                crt(r1, m1, rs[1], ms[1]), m1*ms[1],
                view(rs, 2:n), view(ms, 2:n)
    )
end

function reconstruct_crt(gbs, moduli, Rzz)
    modulo = prod(BigInt.(moduli))
    modelgb = first(gbs)

    reconstructed = zeros(Rzz, length(modelgb))

    for poly_i in 1:length(modelgb)
        builder = MPolyBuildCtx(Rzz)
        for monom in exponent_vectors(modelgb[poly_i])
            mcoeffs = [
                Int(lift(coeff(fs[poly_i], monom)))
                for fs in gbs
            ]
            newcoeff = BigInt(crt(mcoeffs, moduli))
            push_term!(builder, newcoeff, monom)
        end
        reconstructed[poly_i] = finish(builder)
    end

    reconstructed, modulo
end

function reconstruct_modulo(fs, modulo, Rqq)
    ans = zeros(Rqq, length(fs))
    for i in 1:length(fs)
        coerced = map_coefficients(c -> rational_reconstruction(c, modulo), fs[i])
        ans[i] = change_base_ring(base_ring(Rqq), coerced, parent=Rqq)
    end
    ans
end

function run_modular(ps)
    ans = QQ(0)

    for (a, m) in ps
        ans += rational_reconstruction_2(a, m)
        ans -= rational_reconstruction_2(a, m)
    end

    return ans
end
