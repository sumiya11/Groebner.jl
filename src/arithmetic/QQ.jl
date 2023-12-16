# Arithmetic in the rationals.

# All implementations of arithmetic in Q are a subtype of this 
abstract type AbstractArithmeticQQ{AccumType, CoeffType} <:
              AbstractArithmetic{AccumType, CoeffType} end

# Standard arithmetic that uses Base.GMP.MPQ
struct BuiltinArithmeticQQ{AccumType, CoeffType} <:
       AbstractArithmeticQQ{AccumType, CoeffType}
    buf1::BigInt
    buf2::Rational{BigInt}
    function BuiltinArithmeticQQ()
        new{Rational{BigInt}, Rational{BigInt}}(BigInt(0), Rational{BigInt}(0))
    end
end

function select_arithmetic(
    characteristic,
    ::Type{CoeffType},
    hint::Symbol,
    _
) where {CoeffType <: CoeffQQ}
    BuiltinArithmeticQQ()
end
