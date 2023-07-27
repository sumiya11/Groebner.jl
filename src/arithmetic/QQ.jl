# Arithmetic in the rationals.

# All implementations of arithmetic in Q are a subtype of this 
abstract type AbstractArithmeticQQ <: AbstractArithmetic end

# Standard arithmetic that uses Base.GMP.MPQ
struct BuiltinArithmeticQQ <: AbstractArithmeticQQ
    buf1::BigInt
    buf2::Rational{BigInt}
    function BuiltinArithmeticQQ()
        new(BigInt(0), Rational{BigInt}(0))
    end
end

function select_arithmetic(coeffs::Vector{Vector{T}}, ch) where {T <: CoeffQQ}
    BuiltinArithmeticQQ()
end
