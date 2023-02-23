import SIMD: LVec
import SIMD.Intrinsics: llvm_type

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

@inline function baseline_muladd_mod(b::L, x::L, y::L, magic) where {L}
    (b + x*y) % magic
end

p = UInt64(2^31-1)
b = UInt64(2^20 + 1299)
x = UInt64(2^29)
y = UInt64(2^29 + 3)
magic = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(p)
baseline_muladd_mod(b, x, y, magic)

@code_native debuginfo=:none baseline_muladd_mod(b, x, y, magic)

# b + x*y (mod p)
@inline function barrett_muladd_mod(b::U, x::U, y::U, ϕ::L, p::U) where {U, L}
    c = ((x * UInt64(ϕ)) >> (32)) % U
    d = b + x*y - c*p
    min(d - p, d)
end

@inline function barrett_muladd_mod_half(b::U, x::U, y::U, ϕ::L, p::U) where {U, L}
    c = (x * ϕ) >> (32)
    b + x*y - c*p
end

@inline function barrett_muladd_mod_half_v(b, x, y, ϕ, p, s) where {U, L}
    c = (x * ϕ) >> s
    b + x*y - c*p
end

m = 30
@assert m <= 8*sizeof(UInt64)
p = UInt64(Primes.nextprime(2^29))
@assert p <= UInt(2)^(8*sizeof(UInt64)-2)

b = UInt64(2^20 + 1299)
x = UInt64(2^28 + 22)
y = UInt64(2^28 + 3)
ϕ = UInt64(div(UInt64(y) << (sizeof(UInt32)*8), p, RoundUp))
magic = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(p)

baseline_muladd_mod(b, x, y, magic)
barrett_muladd_mod(b, x, y, ϕ, p)
barrett_muladd_mod_half(b, x, y, ϕ, p)

@code_native debuginfo=:none barrett_muladd_mod(b, x, y, ϕ, p)

@code_native debuginfo=:none barrett_muladd_mod_half(b, x, y, ϕ, p)

# X + Y*Z modulo p
function montgomery_2(X, Y, Z, chi, mu, p, m)
    YZ = Y * Z
    b = (YZ * chi) & mu
    c_hi = 1
end

_shuffle_vec(I) = join((string("i32 ", i == :undef ? "undef" : Int32(i::Integer)) for i in I), ", ")

_unpack_hi_bitcast_to(::Type{LVec{2,UInt64}}) = LVec{4,UInt32}
_unpack_hi_bitcast_to(::Type{LVec{4,UInt64}}) = LVec{8,UInt32}
_unpack_hi_shuffle_indices(::Type{LVec{2,UInt64}}) = (1, 4, 3, 5)
_unpack_hi_shuffle_indices(::Type{LVec{4,UInt64}}) = (2, 10, 3, 11, 6, 14, 7, 15)

@generated function _unpack_hi(x::T, y::T) where {T}
    inds = _unpack_hi_shuffle_indices(T)
    shfl = _shuffle_vec(inds)
    U = _unpack_hi_bitcast_to(T)
    M = length(inds)
    llvm = """
        %3 = bitcast $(llvm_type(T)) %0 to $(llvm_type(U))
        %4 = bitcast $(llvm_type(T)) %1 to $(llvm_type(U))
        %5 = shufflevector $(llvm_type(U)) %3, $(llvm_type(U)) %4, <$M x i32> <$shfl>
        %6 = bitcast $(llvm_type(U)) %5 to $(llvm_type(T))
        ret $(llvm_type(T)) %6
        """
    return :(
        $(Expr(:meta, :inline));
        Base.llvmcall($llvm, T, Tuple{T, T}, x, y)
    )
end

@inline function unpack_hi(x::SIMD.Vec{N, T}, y::SIMD.Vec{N, T}) where {N, T}
    SIMD.Vec(_unpack_hi(x.data, y.data))
end

a = SIMD.Vec{2, UInt64}((0x0003fffffffffff9, 0x12300000000fff88))
zero = SIMD.Vec{2, UInt64}((0, 0))

unpack_hi(a, zero)

@code_llvm unpack_hi(a, zero)

