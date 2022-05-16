

# Rational number reconstruction implementation borrowed from CLUE
# and modified a bit to suit the 'Modern Computer Algebra' definitions
# Returns a rational r // h of QQ field in a canonical form such that
#   r // h ≡ a (mod m)
#
# let n = max( λ(a), λ(m) ) , where λ(x) is a number of bits for x
# O(n^2)
function rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}
    bnd = sqrt(float(m) / 2)

    steps = 0

    # @warn "" bnd
    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T

        steps += 1
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])

    # @warn "" r t

    # changed from `<= bnd` to `<= m / bnd`
    # we can speed up this !
    if t <= bnd && gcd(r, t) == 1
        # return (r, t)
        return (true, steps)
    end

    return (false, steps)
end


function maximal_quotient_rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}

    steps = 0

    T = 2^10 * log2(m)

    n, d = 0, 0

    t0, r0 = 0, m
    t1, r1 = 1, a

    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while r1 != 0 && r0 > T
        q = div(r0, r1)

        if q > T
            n, d, T = r1, t1, q
        end

        steps += 1

        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1
    end

    if d == 0 || gcd(n, d) != 1
        return (false, steps)
    end

    if d < 0
        n, d = -n, -d
    end

    return (true, steps)
end
