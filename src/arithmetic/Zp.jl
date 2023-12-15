# Arithmetic in Zp.

# All implementations of arithmetic are a subtype of this
abstract type AbstractArithmetic{AccumType, CoeffType} end

# All implementations of arithmetic in Zp are a subtype of this
abstract type AbstractArithmeticZp{AccumType, CoeffType} <:
              AbstractArithmetic{AccumType, CoeffType} where {AccumType, CoeffType} end

###
# ArithmeticZp

less_than_half(p, ::Type{T}) where {T} = p < (typemax(T) >> ((8 >> 1) * sizeof(T)))

# Modular arithmetic based on builtin classes in `Base.MultiplicativeInverses`
struct ArithmeticZp{AccumType, CoeffType} <: AbstractArithmeticZp{AccumType, CoeffType}
    # magic contains the precomputed multiplicative inverse of the divisor
    magic::UnsignedMultiplicativeInverse{AccumType}

    function ArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        @invariant less_than_half(p, AccumType)
        new{AccumType, CoeffType}(
            UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        )
    end
end

divisor(arithm::ArithmeticZp) = arithm.magic.divisor

mod_p(a::T, arithm::ArithmeticZp{T}) where {T} = a % arithm.magic

inv_mod_p(a::T, arithm::ArithmeticZp{T}) where {T} = invmod(a, divisor(arithm))

###
# SpecializedArithmeticZp

# Same as ArithmeticZp, but stores the fields inline and additionally
# specializes on magic.add
struct SpecializedArithmeticZp{AccumType, CoeffType, Add} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    shift::UInt8
    divisor::AccumType

    function SpecializedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        @invariant less_than_half(p, T)
        uinv = UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        SpecializedArithmeticZp(AccumType, CoeffType, uinv)
    end

    function SpecializedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        uinv::UnsignedMultiplicativeInverse{AccumType}
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        new{AccumType, CoeffType, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
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

# NOTE: the compiler usually fails to simd this for T = UInt32 and T = UInt64.
# However, with T = UInt8 and T = UInt16, it uses the simd mulhi instructions
# such as pmulhuw.
# Our attempts to emulate pmulhuw for T = UInt32 while preserving the
# performance and portability were not successful.
#
# a modulo p (addition specialization)
function mod_p(a::T, mod::SpecializedArithmeticZp{T, C, true}) where {T, C}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
function mod_p(a::T, mod::SpecializedArithmeticZp{T, C, false}) where {T, C}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

inv_mod_p(a::T, arithm::SpecializedArithmeticZp{T}) where {T} = invmod(a, divisor(arithm))

###
# DelayedArithmeticZp

# Same as SpecializedArithmeticZp, but also exploits the case when the
# representation of the modulo has spare leading bits
struct DelayedArithmeticZp{AccumType, CoeffType, Add} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    shift::UInt8
    divisor::AccumType

    function DelayedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        @invariant less_than_half(p, T)
        uinv = UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        DelayedArithmeticZp(AccumType, CoeffType, uinv)
    end

    function DelayedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        uinv::UnsignedMultiplicativeInverse{AccumType}
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        new{AccumType, CoeffType, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
end

divisor(arithm::DelayedArithmeticZp) = arithm.divisor

function n_spare_bits(arithm::DelayedArithmeticZp{T}) where {T <: Unsigned}
    res = leading_zeros(divisor(arithm)) - (8 >> 1) * sizeof(T)
    @invariant res >= 0
    res
end

function n_safe_consecutive_additions(arithm::DelayedArithmeticZp{T}) where {T <: Unsigned}
    T(1) << (2 * n_spare_bits(arithm) - 1)
end

# a modulo p (addition specialization)
function mod_p(a::T, mod::DelayedArithmeticZp{T, true}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
function mod_p(a::T, mod::DelayedArithmeticZp{T, false}) where {T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

inv_mod_p(a::T, arithm::DelayedArithmeticZp{T}) where {T} = invmod(a, divisor(arithm))

###
# CompositeArithmeticZp

# Operates on several integers at once
struct CompositeArithmeticZp{AccumType, CoeffType, TT} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    arithmetics::TT

    function CompositeArithmeticZp(
        ::CompositeInt{2, AccumType},
        ::CompositeInt{2, CoeffType},
        p::CompositeInt{2, CoeffType}
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        a1 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[1])
        a2 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[2])
        new{
            CompositeInt{2, AccumType},
            CompositeInt{2, CoeffType},
            Tuple{typeof(a1), typeof(a2)}
        }((a1, a2))
    end

    function CompositeArithmeticZp(
        ::CompositeInt{4, AccumType},
        ::CompositeInt{4, CoeffType},
        p::CompositeInt{4, CoeffType}
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        a1, a2, a3, a4 = map(SpecializedArithmeticZp, p.data)
        new{
            CompositeInt{4, AccumType},
            CompositeInt{4, CoeffType},
            Tuple{typeof(a1), typeof(a2), typeof(a3), typeof(a4)}
        }((a1, a2, a3, a4))
    end
end

divisor(arithm::CompositeArithmeticZp) = CompositeInt(map(divisor, arithm.arithmetics))

function mod_p(a::CompositeInt{2, T}, arithm::CompositeArithmeticZp) where {T}
    a1 = mod_p(a.data[1], arithm.arithmetics[1])
    a2 = mod_p(a.data[2], arithm.arithmetics[2])
    CompositeInt((a1, a2))
end

function mod_p(a::CompositeInt{4, T}, arithm::CompositeArithmeticZp) where {T}
    a1 = mod_p(a.data[1], arithm.arithmetics[1])
    a2 = mod_p(a.data[2], arithm.arithmetics[2])
    a3 = mod_p(a.data[3], arithm.arithmetics[3])
    a4 = mod_p(a.data[4], arithm.arithmetics[4])
    CompositeInt((a1, a2, a3, a4))
end

function inv_mod_p(a::T, arithm::CompositeArithmeticZp) where {T}
    CompositeInt(invmod.(a.data, divisor(arithm).data))
end

###
# SignedArithmeticZp

struct SignedArithmeticZp{AccumType, CoeffType} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    p::AccumType
    p2::AccumType # = p*p
    magic::Base.MultiplicativeInverses.SignedMultiplicativeInverse{AccumType}

    function SignedArithmeticZp(
        ::Type{AccumType},
        Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffFF, CoeffType <: CoeffFF}
        pa = convert(AccumType, p)
        magic = Base.MultiplicativeInverses.SignedMultiplicativeInverse{AccumType}(pa)
        new{AccumType, CoeffType}(pa, pa * pa, magic)
    end
end

divisor(arithm::SignedArithmeticZp) = arithm.p

function mod_p(a::T, arithm::SignedArithmeticZp{T}) where {T}
    res = a % arithm.magic
    ifelse(res > divisor(arithm), res - divisor(arithm), res)
end

function inv_mod_p(a::T, arithm::SignedArithmeticZp{T}) where {T}
    res = invmod(a, divisor(arithm))
    ifelse(res > divisor(arithm), res - divisor(arithm), res)
end

###
# Selection of arithmetic

"""
    select_arithmetic

Given the `characteristic` of the ground field Z/Zp, the type of the
coefficients `CoeffType`, and the `hint` from the user, returns the most
suitable algorithm for doing arithmetic in Z/Zp.
"""
function select_arithmetic(
    characteristic::Integer,
    ::Type{CoeffType},
    hint::Symbol,
    store_coeffs_tight::Bool
) where {CoeffType <: Union{CoeffFF, CompositeCoeffFF}}
    # NOTE: characteristic is guaranteed to be representable by CoeffType. Maybe
    # change the type of characteristic to CoeffType?
    @assert characteristic < typemax(CoeffType)

    # The type that would act as an accumulator for coefficients of type
    # CoeffType. Usually, this type should be a bit wider than CoeffType
    AccumType = if store_coeffs_tight
        # If the coefficients are stored tightly, 
        # say, using UInt32 with characteristic = 2^31-1
        widen(CoeffType)
    else
        CoeffType
    end

    if hint === :signed
        return SignedArithmeticZp(AccumType, CoeffType, convert(AccumType, characteristic))
    end

    if hint === :delayed
        return DelayedArithmeticZp(AccumType, CoeffType, convert(AccumType, characteristic))
    end

    if hint === :basic
        SpecializedArithmeticZp(AccumType, CoeffType, convert(AccumType, characteristic))
    end

    if hint === :auto
        if CoeffType <: CompositeCoeffFF
            return CompositeArithmeticZp(AccumType, CoeffType, characteristic)
        end

        @assert CoeffType <: Unsigned

        # Use delayed modular arithmetic if the prime is one of the following:
        #   2^16 < prime < 2^29
        #   2^8  < prime < 2^13
        #          prime < 2^5
        # The trade-offs are a bit shifted for large moduli that are > 2^32. It
        # looks like it rarely pays off to use delayed modular arithmetic in
        # such cases
        if 2 < (leading_zeros(characteristic) - (8 >> 1) * sizeof(AccumType))
            return DelayedArithmeticZp(
                AccumType,
                CoeffType,
                convert(AccumType, characteristic)
            )
        else
            return SpecializedArithmeticZp(
                AccumType,
                CoeffType,
                convert(AccumType, characteristic)
            )
        end
    end

    SpecializedArithmeticZp(AccumType, CoeffType, convert(AccumType, characteristic))
end
