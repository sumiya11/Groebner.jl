# Arithmetic in Zp.

# All implementations of arithmetic in Zp are a subtype of it. 
abstract type AbstractArithmeticZp end

# Modular arithmetic implemented based on the Julia builtin classes in
# `Base.MultiplicativeInverses`. Used for primes smaller than 2^32.
struct BuiltinArithmeticZp{T <: Unsigned} <: AbstractArithmeticZp
    # magic contains is a precomputed multiplicative inverse of the divisor
    magic::UnsignedMultiplicativeInverse{T}
    function BuiltinArithmeticZp(p::T) where {T <: Unsigned}
        new{T}(UnsignedMultiplicativeInverse{T}(p))
    end
end

divisor(arithm::BuiltinArithmeticZp) = arithm.magic.divisor

# Same as the built-in implementation, by specializes on the type of the prime
# number being used and stores the fields inline. This implementation is
# preferred for primes up to 128 bit.
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
@inline function mod_x(a::T, mod::SpecializedBuiltinArithmeticZp{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
@inline function mod_x(a::T, mod::SpecializedBuiltinArithmeticZp{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end


function select_arithmetic(coeffs::Vector{Vector{T}}, ch) where {T <: CoeffFF}
    SpecializedBuiltinArithmeticZp(convert(T, ch))
end

# arithmetic over rational numbers
function select_arithmetic(coeffs::Vector{Vector{T}}, ch) where {T <: CoeffQQ}
    ch
end
