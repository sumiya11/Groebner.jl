# Modular arithmetic -- arithmetic in Zp.

# The file provides implementation(s) of modular arithmetic. 
# Modular arithmetic is used in f4/matrix.jl and f4/f4.jl in case
# some of the modular computation techniques are involved.
# It is really important to do modular computation as fast as possible.

# AbstractModularArithmetic. All implementations are a subtype of it. 
abstract type AbstractModularArithmetic end

#------------------------------------------------------------------------------

# Modular arithmetic implemented based on the Julia builtin classes in `Base.MultiplicativeInverses`.
# Used for primes smaller than 2^32.
struct BuiltinModularArithmetic{T <: Unsigned} <: AbstractModularArithmetic
    magic::Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{T}

    function BuiltinModularArithmetic(p::T) where {T <: Unsigned}
        new{T}(Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{T}(p))
    end
end

divisor(mod::BuiltinModularArithmetic) = mod.magic.divisor

# x modulo p
@inline function mod_x(x::T, mod::BuiltinModularArithmetic{T}) where {T}
    x % mod.magic
end

# (a + b*c) modulo p
@inline function mod_muladd(a::T, b::T, c::T, mod::BuiltinModularArithmetic{T}) where {T}
    (a + b * c) % mod.magic
end

# (a*b) modulo p
@inline function mod_mul(a::T, b::T, mod::BuiltinModularArithmetic{T}) where {T}
    (a * b) % mod.magic
end

#------------------------------------------------------------------------------

# Same as the built-in implementation, 
# by specializes (just a bit) on the type of the prime number being used.
# This implementation is preferred for UInt128.
struct SpecializedBuiltinModularArithmetic{T <: Unsigned, Add} <: AbstractModularArithmetic
    multiplier::T
    shift::UInt8
    divisor::T

    function SpecializedBuiltinModularArithmetic(
        uinv::Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{T}
    ) where {T <: Unsigned}
        new{T, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
    function SpecializedBuiltinModularArithmetic(p::T) where {T <: Unsigned}
        uinv = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(p)
        SpecializedBuiltinModularArithmetic(uinv)
    end
end

divisor(mod::SpecializedBuiltinModularArithmetic) = mod.divisor

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

# a modulo p (addition specialization)
@inline function mod_x(a::T, mod::SpecializedBuiltinModularArithmetic{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
@inline function mod_x(a::T, mod::SpecializedBuiltinModularArithmetic{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end
