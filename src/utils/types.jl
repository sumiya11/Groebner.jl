# Common types used throughout the project

# MonomialDegreeOverflow is thrown if there is a risk of monomial degree
# overflow. If we catch a MonomialDegreeOverflow, there is some hope to recover
# the program by restarting with a wider integer type for storing exponents.
struct MonomialDegreeOverflow <: Exception
    msg::String
end

Base.showerror(io::IO, e::MonomialDegreeOverflow) = print(io, e.msg)

###
# All supported coefficient types

# CoeffFF is a supertype of polynomial coefficients modulo a prime
const CoeffFF = Union{Unsigned, Signed}

const CompositeCoeffFF = CompositeInt{N, T} where {N, T <: CoeffFF}

# CoeffQQ is a supertype of polynomial coefficients in the rationals
const CoeffQQ = Union{Rational{BigInt}}

# CoeffZZ is a supertype of polynomial coefficients in the integers.
# NOTE: it is not actually possible to run our implementation of F4 over BigInts
# (or over any other ring that is not a field, in fact). This option should
# probably be deleted
const CoeffZZ = Union{BigInt}

# All supported coefficient types in F4.
const Coeff = Union{CoeffFF, CompositeCoeffFF, CoeffQQ, CoeffZZ}

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
