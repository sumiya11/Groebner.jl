# Arithmetic in the rationals.

# All implementations of arithmetic in QQ are a subtype of it. 
abstract type AbstractArithmeticQQ end

struct BuiltinArithmeticQQ <: AbstractArithmeticQQ
    buf1::BigInt
    buf2::Rational{BigInt}
    function BuiltinArithmeticQQ()
        new(BigInt(0), Rational{BigInt}(0))
    end
end

# arithmetic over rational numbers
function select_arithmetic(coeffs::Vector{Vector{T}}, ch) where {T <: CoeffQQ}
    BuiltinArithmeticQQ()
end
