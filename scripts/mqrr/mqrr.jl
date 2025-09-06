
# Rational reconstruction
function mqrr(x::BigInt, m::BigInt)
    if m <= 0
        throw(ArgumentError("modulus must be positive"))
    end
    if x < 0 || x >= m
        throw(ArgumentError("x must be in [0, m)"))
    end
    T = 2^20*ceil(Int, log(2, m))
    t0, r0 = BigInt(0), m
    t1, r1 = BigInt(1), x
    n, d = 0, 0
    while r1 != 0 && r0 > T
        q = div(r0, r1)
        if q > T
            n, d, T = r1, t1, q
        end
        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1
    end
    (d == 0 || gcd(n, d) != 1) && return false, n, d
    if d < 0
        n, d = -n, -d
    end
    true, n, d
end
