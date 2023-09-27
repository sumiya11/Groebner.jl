# Some common types used throughout the code

# F4 supports several coefficient types:
#   1. Integers modulo a prime, 
#   2. Rational numbers.
# and several exponent-vector implementations:
#   1. Vector, 
#   2. Packed vector,
#   3. Sparse vector.

# MonomialDegreeOverflow is thrown if there is a risk of monomial degree
# overflow. If we catch a MonomialDegreeOverflow, there is some hope to recover
# the program by restarting with a wider int type for exponents.
struct MonomialDegreeOverflow <: Exception
    msg::String
end

Base.showerror(io::IO, e::MonomialDegreeOverflow) = print(io, e.msg)

###
# CoeffFF is a union of all types that a polynomial coefficient modulo a prime
# can have
if isdefined(Groebner, :UInt128)
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
else
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64}
end

# CoeffQQ is a union of all types that a polynomial coefficient in the rationals
# can have
const CoeffQQ = Union{Rational{BigInt}}

# CoeffZZ is a union of all types that a polynomial coefficient in the integers
# can have
# NOTE: it is not actually possible to run F4 over BigInts (or over any other
# ring, in fact). This option is left in compliance with some internals of
# Groebner
const CoeffZZ = Union{BigInt}

# All supported coefficient types in F4.
const Coeff = Union{CoeffFF, CoeffQQ, CoeffZZ}

# Coefficient type used in a single run of classic modular computation
const CoeffModular = UInt64
@assert CoeffModular <: CoeffFF

###
# All supported monomial implementations in F4
const Monom = Union{
    ExponentVector{T} where {T},
    AbstractPackedTuple,
    SparseExponentVector{T, N, I} where {T, N, I}
}
