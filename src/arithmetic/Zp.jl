# Arithmetic in Zp.

# All implementations of arithmetic are a subtype of this
abstract type AbstractArithmetic end

# All implementations of arithmetic in Zp are a subtype of this
abstract type AbstractArithmeticZp <: AbstractArithmetic end

###
# ArithmeticZp

less_than_half(p, ::Type{T}) where {T} = p < (typemax(T) >> ((8 >> 1) * sizeof(T)))

# Modular arithmetic based on builtin classes in `Base.MultiplicativeInverses`
struct ArithmeticZp{T <: Unsigned} <: AbstractArithmeticZp
    # magic contains the precomputed multiplicative inverse of the divisor
    magic::UnsignedMultiplicativeInverse{T}

    function ArithmeticZp(p::T) where {T <: Unsigned}
        @invariant less_than_half(p, T)
        new{T}(UnsignedMultiplicativeInverse{T}(p))
    end
end

divisor(arithm::ArithmeticZp) = arithm.magic.divisor

###
# SpecializedArithmeticZp

# Same as ArithmeticZp, but stores the fields inline and additionally
# specializes on magic.add
struct SpecializedArithmeticZp{T <: Unsigned, Add} <: AbstractArithmeticZp
    multiplier::T
    shift::UInt8
    divisor::T

    function SpecializedArithmeticZp(p::T) where {T <: Unsigned}
        @invariant less_than_half(p, T)
        uinv = UnsignedMultiplicativeInverse(p)
        SpecializedArithmeticZp(uinv)
    end

    function SpecializedArithmeticZp(
        uinv::UnsignedMultiplicativeInverse{T}
    ) where {T <: Unsigned}
        new{T, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
end

divisor(arithm::SpecializedArithmeticZp) = arithm.divisor

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
function mod_p(a::T, mod::SpecializedArithmeticZp{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
function mod_p(a::T, mod::SpecializedArithmeticZp{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

###
# UseSpareBitsArithmeticZp

# Same as SpecializedArithmeticZp, but also exploits the case when the
# representation of the squared modulo in type T has spare leading bits
struct UseSpareBitsArithmeticZp{T <: Unsigned, Add} <: AbstractArithmeticZp
    multiplier::T
    shift::UInt8
    divisor::T

    function UseSpareBitsArithmeticZp(p::T) where {T <: Unsigned}
        @invariant less_than_half(p, T)
        uinv = UnsignedMultiplicativeInverse(p)
        UseSpareBitsArithmeticZp(uinv)
    end

    function UseSpareBitsArithmeticZp(
        uinv::UnsignedMultiplicativeInverse{T}
    ) where {T <: Unsigned}
        new{T, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
end

divisor(arithm::UseSpareBitsArithmeticZp) = arithm.divisor

function n_spare_bits(arithm::UseSpareBitsArithmeticZp{T}) where {T <: Unsigned}
    res = leading_zeros(divisor(arithm)) - (8 >> 1) * sizeof(T)
    @invariant res >= 0
    res
end

function n_safe_consecutive_additions(
    arithm::UseSpareBitsArithmeticZp{T}
) where {T <: Unsigned}
    T(1) << (2 * n_spare_bits(arithm) - 1)
end

# a modulo p (addition specialization)
function mod_p(a::T, mod::UseSpareBitsArithmeticZp{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
function mod_p(a::T, mod::UseSpareBitsArithmeticZp{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

###
# Selection of arithmetic

function select_arithmetic(
    characteristic::Integer,
    ::Type{CoeffType}
) where {CoeffType <: CoeffFF}
    @assert characteristic < typemax(CoeffType)
    # If more than 2 bits are available, that is, characteristic < 2^29
    if 2 < (leading_zeros(characteristic) - (8 >> 1) * sizeof(CoeffType))
        UseSpareBitsArithmeticZp(convert(CoeffType, characteristic))
    else
        SpecializedArithmeticZp(convert(CoeffType, characteristic))
    end
end
