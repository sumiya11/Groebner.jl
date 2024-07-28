# This file is a part of Groebner.jl. License is GNU GPL v2.

# Arithmetic in Zp.

# The arithmetic type is parametrized by the type of the accumulator, AccumType,
# and the type of the coefficient, CoeffType. Generally,
#   a + b*c
# must be representable in AccumType for feasible a,b,c of type CoeffType. One
# example is AccumType = UInt64 and CoeffType = UInt32 with prime moduli.
abstract type AbstractArithmetic{AccumType, CoeffType} end

# All implementations of arithmetic in Z_p are a subtype of this.
#
# Implementations of AbstractArithmeticZp need to implement three functions to
# be usable in F4:
# - divisor(arithmetic)      : returns the prime number p
# - mod_p(x, arithmetic)     : returns x modulo p
# - inv_mod_p(x, arithmetic) : returns x^(-1) modulo p
abstract type AbstractArithmeticZp{AccumType, CoeffType} <:
              AbstractArithmetic{AccumType, CoeffType} end

###
# ArithmeticZp

less_than_half(p, ::Type{T}) where {T} = p < (typemax(T) >> ((8 >> 1) * sizeof(T)))

# Modular arithmetic based on builtin classes in Base.MultiplicativeInverses. 
struct ArithmeticZp{AccumType, CoeffType} <: AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    shift::UInt8
    divisor::AccumType
    add::Bool

    ArithmeticZp(::Type{A}, ::Type{C}, p) where {A, C} = ArithmeticZp(A, C, C(p))

    function ArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant less_than_half(p, AccumType)
        @invariant Primes.isprime(p)
        uinv = UnsignedMultiplicativeInverse{AccumType}(convert(AccumType, p))
        # Further in the code we need the guarantee that the shift is < 64
        @invariant uinv.shift < 8 * sizeof(AccumType)
        new{AccumType, CoeffType}(uinv.multiplier, uinv.shift, uinv.divisor, uinv.add)
    end
end

divisor(arithm::ArithmeticZp) = arithm.divisor

@inline function mod_p(a::T, mod::ArithmeticZp{T}) where {T}
    x = _mul_high(a, mod.multiplier)
    x = ifelse(mod.add, convert(T, convert(T, (convert(T, a - x) >>> UInt8(1))) + x), x)
    unsafe_assume(mod.shift < 8 * sizeof(T))
    a - (x >>> mod.shift) * mod.divisor
end

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
        # Further in the code we need the guarantee that the shift is < 64
        @invariant uinv.shift < 8 * sizeof(AccumType)
        new{AccumType, CoeffType, uinv.add}(uinv.multiplier, uinv.shift, uinv.divisor)
    end
end

divisor(arithm::SpecializedArithmeticZp) = arithm.divisor

# Returns the higher half of the product a*b
@inline function _mul_high(a::T, b::T) where {T <: Union{Signed, Unsigned}}
    ((widen(a) * b) >>> (sizeof(a) * 8)) % T
end

@inline function _mul_high(a::UInt128, b::UInt128)
    shift = sizeof(a) * 4
    mask = typemax(UInt128) >> shift
    a1, a2 = a >>> shift, a & mask
    b1, b2 = b >>> shift, b & mask
    a1b1, a1b2, a2b1, a2b2 = a1 * b1, a1 * b2, a2 * b1, a2 * b2
    carry = ((a1b2 & mask) + (a2b1 & mask) + (a2b2 >>> shift)) >>> shift
    a1b1 + (a1b2 >>> shift) + (a2b1 >>> shift) + carry
end

# NOTE: the compiler sometimes fails to simd this for T ≥ UInt32.
@inline function mod_p(a::T, mod::SpecializedArithmeticZp{T, C, true}) where {T, C}
    x = _mul_high(a, mod.multiplier)
    x = convert(T, convert(T, (convert(T, a - x) >>> UInt8(1))) + x)
    # Bit shifts in Julia check that the shift value is less than the bitsize of
    # the argument type. The assumption here allows us to bypass this check.
    unsafe_assume(mod.shift < 8 * sizeof(T))
    a - (x >>> mod.shift) * mod.divisor
end

@inline function mod_p(a::T, mod::SpecializedArithmeticZp{T, C, false}) where {T, C}
    x = _mul_high(a, mod.multiplier)
    unsafe_assume(mod.shift < 8 * sizeof(T))
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
        @invariant uinv.shift < 8 * sizeof(AccumType)
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
@inline function mod_p(a::A, mod::DelayedArithmeticZp{A, T, true}) where {A, T}
    x = _mul_high(a, mod.multiplier)
    x = convert(A, convert(A, (convert(A, a - x) >>> UInt8(1))) + x)
    unsafe_assume(mod.shift < 8 * sizeof(T))
    a - (x >>> mod.shift) * mod.divisor
end
# a modulo p (no addition specialization)
@inline function mod_p(a::A, mod::DelayedArithmeticZp{A, T, false}) where {A, T}
    x = _mul_high(a, mod.multiplier)
    unsafe_assume(mod.shift < 8 * sizeof(T))
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
        ::Type{CompositeNumber{2, AccumType}},
        ::Type{CompositeNumber{2, CoeffType}},
        p::CompositeNumber{2, CoeffType}
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        a1 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[1])
        a2 = SpecializedArithmeticZp(AccumType, CoeffType, p.data[2])
        new{
            CompositeNumber{2, AccumType},
            CompositeNumber{2, CoeffType},
            Tuple{typeof(a1), typeof(a2)}
        }((a1, a2))
    end

    function CompositeArithmeticZp(
        ::Type{CompositeNumber{N, AccumType}},
        ::Type{CompositeNumber{N, CoeffType}},
        p::CompositeNumber{N, CoeffType}
    ) where {N, AccumType <: CoeffZp, CoeffType <: CoeffZp}
        ai = ntuple(i -> SpecializedArithmeticZp(AccumType, CoeffType, p.data[i]), N)
        new{
            CompositeNumber{N, AccumType},
            CompositeNumber{N, CoeffType},
            Tuple{map(typeof, ai)...}
        }(
            ai
        )
    end
end

divisor(arithm::CompositeArithmeticZp) = CompositeNumber(map(divisor, arithm.arithmetics))

# No chance vectorizing this!
function mod_p(a::CompositeNumber{2, T}, arithm::CompositeArithmeticZp) where {T}
    a1 = mod_p(a.data[1], arithm.arithmetics[1])
    a2 = mod_p(a.data[2], arithm.arithmetics[2])
    CompositeNumber((a1, a2))
end

function mod_p(a::CompositeNumber{4, T}, arithm::CompositeArithmeticZp) where {T}
    a1 = mod_p(a.data[1], arithm.arithmetics[1])
    a2 = mod_p(a.data[2], arithm.arithmetics[2])
    a3 = mod_p(a.data[3], arithm.arithmetics[3])
    a4 = mod_p(a.data[4], arithm.arithmetics[4])
    CompositeNumber((a1, a2, a3, a4))
end

function inv_mod_p(a::T, arithm::CompositeArithmeticZp) where {T}
    CompositeNumber(invmod.(a.data, divisor(arithm).data))
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
        @invariant magic.shift < 8 * sizeof(AccumType)
        new{AccumType, CoeffType}(pa, pa * pa, magic.multiplier, magic.addmul, magic.shift)
    end
end

divisor(arithm::SignedArithmeticZp) = arithm.p

@inline function mod_p(a::T, mod::SignedArithmeticZp{T}) where {T}
    x = _mul_high(a, mod.multiplier)
    x += a * mod.addmul
    unsafe_assume(mod.shift < 8 * sizeof(T))
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
    ps::CompositeNumber{N, T}
    p2s::CompositeNumber{N, T}

    multipliers::CompositeNumber{N, T}
    addmuls::CompositeNumber{N, Int8}
    shifts::CompositeNumber{N, UInt8}

    function SignedCompositeArithmeticZp(
        ::Type{CompositeNumber{N, AT}},
        ::Type{CompositeNumber{N, CT}},
        ps::CompositeNumber{N, CT}
    ) where {N, AT <: CoeffZp, CT <: CoeffZp}
        arithms = map(pj -> SignedArithmeticZp(AT, CT, pj), ps.data)
        multipliers = map(a -> a.multiplier, arithms)
        addmuls = map(a -> a.addmul, arithms)
        shifts = map(a -> a.shift, arithms)
        @invariant all(x -> x < 8 * sizeof(AT), shifts)
        p2s = CompositeNumber{N, AT}(ps) * CompositeNumber{N, AT}(ps)
        new{CompositeNumber{N, AT}, CompositeNumber{N, CT}, AT, N}(
            CompositeNumber{N, AT}(ps),
            p2s,
            CompositeNumber(multipliers),
            CompositeNumber(addmuls),
            CompositeNumber(shifts)
        )
    end
end

divisor(arithm::SignedCompositeArithmeticZp) = arithm.ps

@inline function mod_p(
    a::CompositeNumber{N, T},
    arithm::SignedCompositeArithmeticZp{CompositeNumber{N, T}, U, W, N}
) where {N, T, U, W}
    x = _mul_high.(a.data, arithm.multipliers.data)
    x = x .+ a.data .* arithm.addmuls.data
    unsafe_assume(all(arithm.shifts.data .< (8 * sizeof(T))))
    d = (signbit.(x) .+ (x .>> arithm.shifts.data)) .% T
    res = a.data .- d .* arithm.ps.data
    res = ifelse.(res .>= T(0), res, res .+ arithm.ps.data)
    CompositeNumber(res)
end

function inv_mod_p(a::T, arithm::SignedCompositeArithmeticZp{T}) where {T}
    CompositeNumber(invmod.(a.data, arithm.ps.data))
end

###
# FloatingPointArithmeticZp

struct FloatingPointArithmeticZp{AccumType, CoeffType} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::AccumType
    divisor::AccumType

    function FloatingPointArithmeticZp(
        ::Type{AccumType},
        ::Type{CoeffType},
        p::CoeffType
    ) where {AccumType <: CoeffZp, CoeffType <: CoeffZp}
        @invariant AccumType === Float64
        @invariant 0 < p < 2^25  # < 52 / 2 - 1
        @invariant Primes.isprime(Int(p))
        multiplier = 1 / p
        new{AccumType, CoeffType}(multiplier, p)
    end
end

divisor(arithm::FloatingPointArithmeticZp) = arithm.divisor

@inline function mod_p(a::T, mod::FloatingPointArithmeticZp{T, C}) where {T, C}
    @fastmath begin
        b = a * mod.multiplier
        c = floor(b)
        a - mod.divisor * c # may be fused
    end
end

@inline function fma_mod_p(
    a1::T,
    a2::T,
    a3::T,
    mod::FloatingPointArithmeticZp{T, C}
) where {T, C}
    @fastmath begin
        b = fma(a1, a2, a3)
        mod_p(b, mod)
    end
end

inv_mod_p(a::T, arithm::FloatingPointArithmeticZp{T}) where {T} =
    T(invmod(Int(a), Int(divisor(arithm))))

###
# FloatingPointCompositeArithmeticZp

struct FloatingPointCompositeArithmeticZp{AccumType, CoeffType, T, N} <:
       AbstractArithmeticZp{AccumType, CoeffType}
    multiplier::CompositeNumber{N, T}
    divisor::CompositeNumber{N, T}

    function FloatingPointCompositeArithmeticZp(
        ::Type{CompositeNumber{N, AT}},
        ::Type{CompositeNumber{N, CT}},
        ps::CompositeNumber{N, CT}
    ) where {N, AT <: CoeffZp, CT <: CoeffZp}
        @invariant AT === Float64
        @invariant all(0 .< ps.data .< 2^25)  # < 52 / 2 - 1
        @invariant all(Primes.isprime.(Int.(ps.data)))
        multiplier = inv(ps)
        new{CompositeNumber{N, AT}, CompositeNumber{N, CT}, AT, N}(multiplier, ps)
    end
end

divisor(arithm::FloatingPointCompositeArithmeticZp) = arithm.divisor

@inline function mod_p(a::T, mod::FloatingPointCompositeArithmeticZp{T, C}) where {T, C}
    @fastmath begin
        b = a * mod.multiplier
        c = T(floor.(b.data))
        a - mod.divisor * c # may be fused
    end
end

inv_mod_p(a::T, arithm::FloatingPointCompositeArithmeticZp{T}) where {T} =
    T(invmod.(Int.(a.data), Int.(divisor(arithm).data)))

###
# Selection of arithmetic

# Returns the most suitable algorithm for doing arithmetic in the ground field.
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
    @invariant characteristic <= typemax(CoeffType)
    @invariant isinteger(characteristic)

    # The type that would act as an accumulator. Usually, this type must be
    # wider than CoeffType
    AccumType = if using_wide_type_for_coeffs
        CoeffType
    else
        # If the coefficients are stored tightly, say, using UInt32 paired with
        # characteristic = 2^31-1, then we need a bigger accumulator
        widen(CoeffType)
    end

    if CoeffType <: CompositeCoeffZp
        if hint === :signed || CoeffType <: CompositeNumber{N, T} where {N, T <: Signed}
            return SignedCompositeArithmeticZp(
                AccumType,
                CoeffType,
                CoeffType(characteristic)
            )
        elseif hint === :floating
            return FloatingPointCompositeArithmeticZp(
                AccumType,
                CoeffType,
                CoeffType(characteristic)
            )
        else
            return CompositeArithmeticZp(AccumType, CoeffType, CoeffType(characteristic))
        end
    end

    if hint === :floating
        return FloatingPointArithmeticZp(AccumType, CoeffType, CoeffType(characteristic))
    end

    if hint === :signed
        return SignedArithmeticZp(AccumType, CoeffType, characteristic)
    end

    if hint === :delayed
        if iszero(leading_zeros(characteristic) - (8 >> 1) * sizeof(AccumType))
            @log :warn "Cannot use $hint arithmetic with characteristic $characteristic"
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
        #           prime < 2^4
        # The trade-offs are a bit shifted for large moduli that are > 2^32. It
        # looks like it rarely pays off to use delayed modular arithmetic in
        # such cases
        if !(AccumType === UInt128)
            if !using_wide_type_for_coeffs &&
               leading_zeros(CoeffType(characteristic)) > 4 ||
               using_wide_type_for_coeffs &&
               ((8 * sizeof(CoeffType)) >> 1) -
               (8 * sizeof(CoeffType) - leading_zeros(CoeffType(characteristic))) > 4
                return DelayedArithmeticZp(
                    AccumType,
                    CoeffType,
                    convert(AccumType, characteristic)
                )
            end
        end
    end

    SpecializedArithmeticZp(AccumType, CoeffType, characteristic)
end
