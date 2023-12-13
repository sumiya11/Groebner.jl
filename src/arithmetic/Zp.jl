# Arithmetic in Zp.

# All implementations of arithmetic are a subtype of this
abstract type AbstractArithmetic end

# All implementations of arithmetic in Zp are a subtype of this
abstract type AbstractArithmeticZp <: AbstractArithmetic end

# Modular arithmetic based on the Julia builtin classes in
# `Base.MultiplicativeInverses`. Used for primes smaller than 2^32.
struct BuiltinArithmeticZp{T <: Unsigned} <: AbstractArithmeticZp
    # magic contains is a precomputed multiplicative inverse of the divisor
    magic::UnsignedMultiplicativeInverse{T}
    function BuiltinArithmeticZp(p::T) where {T <: Unsigned}
        new{T}(UnsignedMultiplicativeInverse{T}(p))
    end
end

divisor(arithm::BuiltinArithmeticZp) = arithm.magic.divisor

function n_reserved_bits(arithm::BuiltinArithmeticZp{T}) where {T <: Unsigned}
    res = 1 + leading_zeros(divisor(arithm)) - (8 >> 1) * sizeof(T)
    @invariant res >= 0
    res
end

function skip(arithm::BuiltinArithmeticZp{T}) where {T <: Unsigned}
    T(1) << (2 * n_reserved_bits(arithm) - 2)
end

# Same as the built-in one, by specializes on the type of the prime number and
# stores the fields inline. This implementation is preferred for primes up to 64
# bits.
struct SpecializedBuiltinArithmeticZp{T <: Unsigned, Add} <: AbstractArithmeticZp
    multiplier::T
    shift::UInt8
    divisor::T
    function SpecializedBuiltinArithmeticZp(
        uinv::UnsignedMultiplicativeInverse{T}
    ) where {T <: Unsigned}
        new{T, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
    function SpecializedBuiltinArithmeticZp(p::T) where {T <: Unsigned}
        uinv = UnsignedMultiplicativeInverse(p)
        SpecializedBuiltinArithmeticZp(uinv)
    end
end

divisor(arithm::SpecializedBuiltinArithmeticZp) = arithm.divisor

function n_reserved_bits(arithm::SpecializedBuiltinArithmeticZp{T}) where {T <: Unsigned}
    # +1 thanks to unsigned representation
    res = 1 + leading_zeros(divisor(arithm)) - (8 >> 1) * sizeof(T)
    @invariant res >= 0
    res
end

function skip(arithm::SpecializedBuiltinArithmeticZp{T}) where {T <: Unsigned}
    T(1) << (2 * n_reserved_bits(arithm) - 2)
end

# Returns the higher half of the product a*b
function _mul_high(a::T, b::T) where {T <: Union{Signed, Unsigned}}
    ((widen(a) * b) >>> (sizeof(a) * 8)) % T
end

function _mul_high(a::UInt128, b::UInt128)
    shift = sizeof(a) * 4
    mask = typemax(UInt128) >> shift
    a1, a2 = a >>> shift, a & mask
    b1, b2 = b >>> shift, b & mask
    a1b1, a1b2, a2b1, a2b2 = a1 * b1, a1 * b2, a2 * b1, a2 * b2
    carry = ((a1b2 & mask) + (a2b1 & mask) + (a2b2 >>> shift)) >>> shift
    a1b1 + (a1b2 >>> shift) + (a2b1 >>> shift) + carry
end

# TODO: move to linear algebra
# a modulo p (addition specialization)
@inline function mod_p(a::T, mod::SpecializedBuiltinArithmeticZp{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
@inline function mod_p(a::T, mod::SpecializedBuiltinArithmeticZp{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

###
# Selection of arithmetic

function select_arithmetic(characteristic, ::Type{CoeffType}) where {CoeffType <: CoeffFF}
    SpecializedBuiltinArithmeticZp(convert(CoeffType, characteristic))
end
