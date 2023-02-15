# Here we define and describe some types used throughout the code.

# F4 supports several coefficient types:
#   1. Integers modulo prime, 
#   2. Integers, 
#   3. Rational numbers
# and several exponent-vector implementations:
#   i. Vector, 
#   ii. Packed vector,
#   iii. Sparse Packed vector(currently not used)

# Integers modulo prime are a subtype of CoeffFF.
# They are implemented as word-sized unsigned integers.

# Among those, UInt64 is the default type used in modular computations.
#=
    The most suitable concrete type is chosen based
    on the size of the field characteristic:

    2^1  <= char < 2^4   ==>   UInt8
    2^4  <= char < 2^8   ==>   UInt16
    2^8  <= char < 2^16  ==>   UInt32
    2^16 <= char < 2^32  ==>   UInt64
    2^32 <= char < 2^64  ==>   UInt128
=#
# Thus, we support prime modulo up to 2^64-1.

# Finite field characteristic will be of the same type
# as the coefficients, hence, a subtype of CoeffFF.

# CoeffFF -- all possible types for finite field computation.
if isdefined(Groebner, :UInt128)
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
else
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64}
end

# CoeffFF -- all possible types for integer ring computation.
const CoeffZZ = Union{BigInt}

# CoeffQQ -- all possible types for rationals computation.
const CoeffQQ = Union{Rational{BigInt}}

# Polynomial coefficient from K[t1,...,tm] 
# (currently not used, to be extended with Gleb!)
const CoeffParam = Union{AbstractAlgebra.PolyElem, AbstractAlgebra.MPolyElem}

# All supported coefficient types in F4
const Coeff = Union{CoeffFF, CoeffZZ, CoeffQQ, CoeffParam}

# Coefficient type used in a single run of modular computation.
# Must be a subtype of CoeffFF
const CoeffModular = UInt64
@assert CoeffModular <: CoeffFF

#------------------------------------------------------------------------------

# Packed monomial exponent vector type;
# Currently, this is implemented as a sequence of word-sized integers.
# which pack a certain number of exponents together.
# In 99% of cases, this representation is used
const PackedExponentVector = Union{AbstractPackedPair}

# Sparse monomial exponent vector type;
# used for computations in very sparse problems
# (currently not used at all)
const SparseExponentVector = Union{Vector}

# Standard monomial exponent vector type.
# Used only as a last resort (in case packed representation fails)
const ExponentVector = Union{PowerVector}

# All supported monomial implementations in F4
const Monom = Union{ExponentVector, SparseExponentVector, PackedExponentVector}

# the type of entry in the monom
powertype(m::Monom) = eltype(typeof(m))
powertype(::Type{M}) where {M<:Monom} = eltype(M)

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

# Index of a monomial in the hashtable
# (see f4/hashtable.jl)
const MonomIdx = Int32

# Hash of a monomial in the hashtable
# !!! Changing a type with another size will cause errors
# (see f4/hashtable.jl)
const MonomHash = UInt32

# Division mask of a monomial
# (see f4/hashtable.jl) 
const DivisionMask = UInt32

const ColumnIdx = Int32

#------------------------------------------------------------------------------

# f4 may fail in some cases and throw a RecoverableException.
# If we catch a RecoverableException, there is a hope to recover the program.
#
# Currently, RecoverableException can be caused by one of the following:
# - overflow in monomial operations (monoms/packedpairs.jl),
# - bad choce of a prime during rational computation (gb/groebner.jl),
# - fail of randomized sparse linear algebra (f4/matrix.jl),
#
struct RecoverableException <: Exception
    msg::String
end

Base.showerror(io::IO, e::RecoverableException) = print(io, e.msg)
