# Arithmetic in Zp.

# All implementations of arithmetic are a subtype of this.
#
# The arithmetic is parametrized by the type of the accumulator and the type of
# the coefficient. Generally, AccumType should be picked so that the result of
#   a + b*c
# is representable exactly in AccumType for feasible  a,b,c of type CoeffType. 
# One common example is AccumType = UInt64 and CoeffType = UInt32
abstract type AbstractArithmetic{AccumType, CoeffType} end

# All implementations of arithmetic in Zp are a subtype of this.
#
# Implementations of AbstractArithmeticZp need to implement three functions to
# be usable in F4:
# - divisor(arithmetic)      : returns the prime number p that defines Z/pZ
# - mod_p(x, arithmetic)     : returns x modulo p
# - inv_mod_p(x, arithmetic) : returns x^(-1) modulo p
abstract type AbstractArithmeticZp{AccumType, CoeffType} <:
              AbstractArithmetic{AccumType, CoeffType} end

###
# ArithmeticZp

less_than_half(p, ::Type{T}) where {T} = p < (typemax(T) >> ((8 >> 1) * sizeof(T)))

# Modular arithmetic based on builtin classes in `Base.MultiplicativeInverses`. 
struct ArithmeticZp{AccumType, CoeffType} <: AbstractArithmeticZp{AccumType, CoeffType}
    # magic contains the precomputed multiplicative inverse of the divisor
    magic::UnsignedMultiplicativeInverse{AccumType}

    ArithmeticZp(::Type{A}, ::Type{C}, p) where {A, C} = ArithmeticZp(A, C, C(p))

    function ArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant less_than_half(p, AccumType)
        @invariant Primes.isprime(p)
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
# NOTE: can specialize even further if we only consider the primes such that
# shift = 0. However, the number of such primes is perhaps < 100 in the range
# 1..2^64 
struct SpecializedArithmeticZp{AccumType, CoeffType, Add} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    shift::UInt8
    divisor::AccumType

    SpecializedArithmeticZp(::Type{A}, ::Type{C}, p) where {A, C} =
        SpecializedArithmeticZp(A, C, C(p))

    function SpecializedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant less_than_half(p, AccumType)
        @invariant Primes.isprime(p)
        uinv = UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        SpecializedArithmeticZp(AccumType, CoeffType, uinv)
    end

    function SpecializedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        uinv::UnsignedMultiplicativeInverse{AccumType}
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
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

# NOTE: the compiler usually fails to simd this operation for T = UInt32 and T =
# UInt64. However, with T = UInt8 and T = UInt16, it uses the simd mulhi
# instructions such as pmulhuw.
# Our attempts to emulate pmulhuw for T = UInt32 while preserving the
# performance and portability were not very promising.
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

# Exploits the case when the representation of the prime has spare leading bits
# and thus delays reduction modulo a prime
struct DelayedArithmeticZp{AccumType, CoeffType, Add} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    shift::UInt8
    divisor::AccumType

    DelayedArithmeticZp(::Type{A}, ::Type{C}, p) where {A, C} =
        DelayedArithmeticZp(A, C, C(p))

    function DelayedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant less_than_half(p, AccumType)
        @invariant leading_zeros(p) > 0
        @invariant Primes.isprime(p)
        uinv = UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        DelayedArithmeticZp(AccumType, CoeffType, uinv)
    end

    function DelayedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        uinv::UnsignedMultiplicativeInverse{AccumType}
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
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
function mod_p(a::A, mod::DelayedArithmeticZp{A, T, true}) where {A, T}
    x = _mul_high(a, mod.multiplier)
    x = convert(A, convert(A, (convert(A, a - x) >>> 1)) + x)
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
function mod_p(a::A, mod::DelayedArithmeticZp{A, T, false}) where {A, T}
    x = _mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

inv_mod_p(a::T, arithm::DelayedArithmeticZp{T}) where {T} = invmod(a, divisor(arithm))

###
# CompositeArithmeticZp

# Can operate on several integers at once
struct CompositeArithmeticZp{AccumType, CoeffType, TT} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    arithmetics::TT

    function CompositeArithmeticZp(
        ::Type{CompositeInt{2, AccumType}},
        ::Type{CompositeInt{2, CoeffType}},
        p::CompositeInt{2, CoeffType}
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        a1 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[1])
        a2 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[2])
        new{
            CompositeInt{2, AccumType},
            CompositeInt{2, CoeffType},
            Tuple{typeof(a1), typeof(a2)}
        }((a1, a2))
    end

    function CompositeArithmeticZp(
        ::Type{CompositeInt{N, AccumType}},
        ::Type{CompositeInt{N, CoeffType}},
        p::CompositeInt{N, CoeffType}
    ) where {N, AccumType <: CoeffZp, CoeffType <: CoeffZp}
        ai = ntuple(i -> SpecializedArithmeticZp(AccumType, CoeffType, p.data[i]), N)
        new{
            CompositeInt{N, AccumType},
            CompositeInt{N, CoeffType},
            Tuple{map(typeof, ai)...}
        }(
            ai
        )
    end
end

divisor(arithm::CompositeArithmeticZp) = CompositeInt(map(divisor, arithm.arithmetics))

# No chance vectorizing this!
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

# Uses signed types and exploits the trick with adding p^2 to negative values
struct SignedArithmeticZp{AccumType, CoeffType} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    p::AccumType
    p2::AccumType # = p*p

    multiplier::AccumType
    addmul::Int8
    shift::UInt8

    SignedArithmeticZp(::Type{A}, ::Type{C}, p) where {A, C} =
        SignedArithmeticZp(A, C, C(p))

    function SignedArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant Primes.isprime(p)
        pa = convert(AccumType, p)
        magic = Base.MultiplicativeInverses.SignedMultiplicativeInverse(pa)
        new{AccumType, CoeffType}(pa, pa * pa, magic.multiplier, magic.addmul, magic.shift)
    end
end

divisor(arithm::SignedArithmeticZp) = arithm.p

function mod_p(a::T, mod::SignedArithmeticZp{T}) where {T}
    x = _mul_high(a, mod.multiplier)
    x += a * mod.addmul
    d = (signbit(x) + (x >> mod.shift)) % T
    res = a - d * mod.p
    ifelse(res >= zero(T), res, res + mod.p)
end

function inv_mod_p(a::T, arithm::SignedArithmeticZp{T}) where {T}
    invmod(a, divisor(arithm))
end

###
# SignedCompositeArithmeticZp

# Operates on several integers at once, uses signed representation
struct SignedCompositeArithmeticZp{AccumType, CoeffType, T, N} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    ps::CompositeInt{N, T}
    p2s::CompositeInt{N, T}

    multipliers::CompositeInt{N, T}
    addmuls::CompositeInt{N, Int8}
    shifts::CompositeInt{N, UInt8}

    function SignedCompositeArithmeticZp(
        ::Type{CompositeInt{N, AT}},
        ::Type{CompositeInt{N, CT}},
        ps::CompositeInt{N, CT}
    ) where {N, AT <: CoeffZp, CT <: CoeffZp}
        arithms = map(pj -> SignedArithmeticZp(AT, CT, pj), ps.data)
        multipliers = map(a -> a.multiplier, arithms)
        addmuls = map(a -> a.addmul, arithms)
        shifts = map(a -> a.shift, arithms)
        p2s = CompositeInt{N, AT}(ps) * CompositeInt{N, AT}(ps)
        new{CompositeInt{N, AT}, CompositeInt{N, CT}, AT, N}(
            CompositeInt{N, AT}(ps),
            p2s,
            CompositeInt(multipliers),
            CompositeInt(addmuls),
            CompositeInt(shifts)
        )
    end
end

divisor(arithm::SignedCompositeArithmeticZp) = arithm.ps

function mod_p(
    a::CompositeInt{N, T},
    arithm::SignedCompositeArithmeticZp{CompositeInt{N, T}, U, W, N}
) where {N, T, U, W}
    x = _mul_high.(a.data, arithm.multipliers.data)
    x = x .+ a.data .* arithm.addmuls.data
    d = (signbit.(x) .+ (x .>> arithm.shifts.data)) .% T
    res = a.data .- d .* arithm.ps.data
    res = ifelse.(res .>= T(0), res, res .+ arithm.ps.data)
    CompositeInt(ifelse.(res .> arithm.ps.data, res .- arithm.ps.data, res))
end

function inv_mod_p(a::T, arithm::SignedCompositeArithmeticZp{T}) where {T}
    CompositeInt(invmod.(a.data, arithm.ps.data))
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
    ::Type{CoeffType},
    characteristic::CharType,
    hint::Symbol,
    using_wide_type_for_coeffs::Bool
) where {
    CoeffType <: Union{CoeffZp, CompositeCoeffZp},
    CharType <: Union{CoeffZp, CompositeCoeffZp}
}
    # We guarantee that the characteristic is representable by CoeffType. 
    # Maybe change the type of characteristic to CoeffType?
    @assert characteristic <= typemax(CoeffType)

    # The type that would act as an accumulator. Usually, this type should be a
    # bit wider than CoeffType
    AccumType = if using_wide_type_for_coeffs
        CoeffType
    else
        # If the coefficients are stored tightly, say, using UInt32 paired with
        # characteristic = 2^31-1, then we need a bigger accumulator
        widen(CoeffType)
    end

    if CoeffType <: CompositeCoeffZp
        if hint === :signed || CoeffType <: CompositeInt{N, T} where {N, T <: Signed}
            return SignedCompositeArithmeticZp(
                AccumType,
                CoeffType,
                CoeffType(characteristic)
            )
        else
            return CompositeArithmeticZp(AccumType, CoeffType, CoeffType(characteristic))
        end
    end

    if hint === :signed
        return SignedArithmeticZp(AccumType, CoeffType, characteristic)
    end

    if hint === :delayed
        if iszero(leading_zeros(characteristic) - (8 >> 1) * sizeof(AccumType))
            @log level = 1_000 "Cannot use $hint arithmetic with characteristic $characteristic"
            @assert false
        end
        return DelayedArithmeticZp(AccumType, CoeffType, characteristic)
    end

    if hint === :basic
        SpecializedArithmeticZp(AccumType, CoeffType, characteristic)
    end

    if hint === :auto
        @assert CoeffType <: Unsigned

        # Use delayed modular arithmetic if the prime is one of the following:
        #   2^16 <= prime < 2^28
        #   2^8  <= prime < 2^12
        #          prime < 2^4
        # The trade-offs are a bit shifted for large moduli that are > 2^32. It
        # looks like it rarely pays off to use delayed modular arithmetic in
        # such cases
        if !(AccumType === UInt128) &&
           4 < (leading_zeros(CoeffType(characteristic)) - (8 >> 1) * sizeof(AccumType))
            return DelayedArithmeticZp(
                AccumType,
                CoeffType,
                convert(AccumType, characteristic)
            )
        else
            return SpecializedArithmeticZp(AccumType, CoeffType, characteristic)
        end
    end

    SpecializedArithmeticZp(AccumType, CoeffType, characteristic)
end
