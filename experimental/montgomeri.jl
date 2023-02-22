# Computing a modulo p.

# a < 2^m p

# 0 < rho < p
# 0 < chi < 2^m
# mu = 2^m - 1
#
# rho 2^m - chi p = 1

# Say,
m = 31
mu = UInt32(2)^m - UInt32(1)
p = UInt32(2^31 - 1)

# Take rho = 1
rho = UInt32(1)
# Then,
chi = div(rho * UInt64(2)^m - 1, p)
# _, rho, chi = gcdx(UInt64(2^m), UInt64(p))

# rho = rho % UInt32
# chi = chi % UInt32

function montgomery_1(a::T, chi, mu, p::T, m) where {T}
    b = (a * chi) & mu
    c = a + widen(b) * p
    d = (c >>> m) % T
    d >= p ? d - p : d
end

a = UInt32(rand(p:2p))
montgomery_1(a, chi, mu, p, m)

a % p

@inline function baseline_mul_mod(x::L, y::L, magic) where {L}
    (x*y) % magic
end

p = UInt64(2^31-1)
x = UInt64(2^29)
y = UInt64(2^29 + 3)
magic = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(p)
baseline_mul_mod(x, y, magic)

function barrett_mul_mod(x::U, y::U, ϕ::L, p::U) where {U, L}
    c = ((x * UInt128(ϕ)) >> (sizeof(U)*8)) % U
    d = x*y - c*p
    min(d - p, d)
end

function barrett_mul_mod_half(x::U, y::U, ϕ::L, p::U) where {U, L}
    c = ((x * UInt128(ϕ)) >> (sizeof(U)*8)) % U
    x*y - c*p
end

p = UInt64(2^31-1)
x = UInt64(2^29)
y = UInt64(2^29 + 3)
ϕ = UInt64(div(UInt128(y) << (sizeof(UInt64)*8), p, RoundUp))
barrett_mul_mod_half(x, y, ϕ, p)

# X + Y*Z modulo p
function montgomery_2(X, Y, Z, chi, mu, p, m)
    YZ = Y * Z
    b = (YZ * chi) & mu
    c_hi = 1
end

a = Vec{4, UInt32}((1, 2, 3, 4))
b = Vec{4, UInt32}((0, 5, 22, 9))
a + b

function pup(x::T, y::T) where {T}
    Base.llvmcall("""
        %3 = add <4 x i32> %1, %0
        ret <4 x i32> %3
        """, T, Tuple{T, T}, x, y)
end

pup(a.data, b.data)
