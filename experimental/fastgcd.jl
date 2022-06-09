
function absdiff(x::T,y::T) where {T<:Unsigned}
    d = max(x,y) - min(x,y)
    d, d
end
function absdiff(x::T,y::T) where {T<:Signed}
    d = x - y
    abs(d), d
end
# binary GCD (aka Stein's) algorithm
# about 1.7x (2.1x) faster for random Int64s (Int128s)
# Unfortunately, we need to manually annotate this as `@assume_effects :terminates_locally` to work around #41694.
# Since this is used in the Rational constructor, constant folding is something we do care about here.
function _gcd(ain::T, bin::T) where T<:Integer
    zb = trailing_zeros(bin)
    za = trailing_zeros(ain)
    a = abs(ain)
    b = abs(bin >> zb)
    k = min(za, zb)
    while a != 0
        a >>= za
        absd, diff = absdiff(a, b)
        za = trailing_zeros(diff)
        b = min(a, b)
        a = absd
    end
    r = b << k
    return r % T
end
