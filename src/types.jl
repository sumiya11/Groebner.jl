# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# All supported monomial implementations in F4

const Monom = Union{ExponentVector{T} where {T}, AbstractPackedTuple}

###
# All supported coefficient types

# CoeffZp is a supertype of polynomial coefficients modulo a prime
const CoeffZp = Union{AbstractFloat, Unsigned, Signed}

const CompositeCoeffZp = CompositeNumber{N, T} where {N, T <: CoeffZp}

# CoeffQQ is a supertype of polynomial coefficients in the rationals
const CoeffQQ = Union{Rational{BigInt}}

# CoeffZZ is a supertype of polynomial coefficients in the integers. Note that
# it is not actually possible to run our F4 over BigInts (or over any other ring
# that is not a field, in fact).
const CoeffZZ = Union{BigInt}

# All supported coefficient types in F4.
const Coeff = Union{CoeffZp, CompositeCoeffZp, CoeffQQ, CoeffZZ, CoeffGeneric}

# Coefficient type used in a single run of classic modular computation
const CoeffModular = UInt64
@assert CoeffModular <: CoeffZp
