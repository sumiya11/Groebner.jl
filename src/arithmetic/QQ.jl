# This file is a part of Groebner.jl. License is GNU GPL v2.

# Arithmetic in the rationals.

# All implementations of arithmetic in the rationals are a subtype of this 
abstract type AbstractArithmeticQQ{AccumType, CoeffType} <: AbstractArithmetic{AccumType, CoeffType} end

# Standard arithmetic that uses Base.GMP.MPQ
struct ArithmeticQQ{AccumType, CoeffType} <: AbstractArithmeticQQ{AccumType, CoeffType}
    buf1::BigInt
    buf2::Rational{BigInt}
    function ArithmeticQQ()
        new{Rational{BigInt}, Rational{BigInt}}(BigInt(0), Rational{BigInt}(0))
    end
end

function select_arithmetic(
    ::Type{CoeffType1},
    characteristic::CoeffType2,
    _,
    _
) where {CoeffType1 <: CoeffQQ, CoeffType2 <: Coeff}
    @invariant iszero(characteristic)
    ArithmeticQQ()
end
