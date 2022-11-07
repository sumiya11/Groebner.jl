
# F4 supports several coefficient types
# (Integers modulo prime, Integers, Rational numbers)
# and several exponent-vector implementations
# (Vector, Packed vector, Sparse Packed vector).
# These are described below in the file.

# Polynomial coefficient from finite field type.
# Among those, UInt64 is the default type used in modular computations
#=
    The most suitable type is chosen based on the size of the
    field characteristic during combination:

    2^1  <= char < 2^4   ==>   UInt8
    2^4  <= char < 2^8   ==>   UInt16
    2^8  <= char < 2^16  ==>   UInt32
    2^16 <= char < 2^32  ==>   UInt64
    2^32 <= char < 2^64  ==>   UInt128
=#

# Finite field characteristic will be of the same type
# as the coefficients, hence, a subtype of CoeffFF
if isdefined(Groebner, :UInt128)
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
else
    const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64}
end

# Polynomial coefficient from integer ring type
const CoeffZZ = Union{BigInt}

# Polynomial coefficient from rational field type
const CoeffQQ = Union{Rational{BigInt}}

# Polynomial coefficient from K[t1,...,tm] 
# to be done with Gleb !
const CoeffParam = Union{AbstractAlgebra.PolyElem, AbstractAlgebra.MPolyElem}

# All supported coefficient types in F4
const Coeff = Union{CoeffFF, CoeffZZ, CoeffQQ, CoeffParam}

#------------------------------------------------------------------------------

# Coefficient type used in a single run of modular computation.
# Must be a subtype of CoeffFF
const CoeffModular = UInt64

@assert CoeffModular <: CoeffFF

#------------------------------------------------------------------------------

# There are several possible exponent vector implementations
# The default one is Packed implementation, which may vary
# based on some heuristics

# Sparse monomial exponent vector type;
# used for computations in very sparse problems
const SparseExponentVector = Union{Vector}

# Packed monomial exponent vector type;
# In 99% of cases, this representation is used
const PackedExponentVector = Union{AbstractPackedPair}

# Standard monomial exponent vector type;
# works for total degrees up to 2^64-1;
# Used as a last resort (in case Packed representation fails)
# (which is the maximum degree possible in Julia)
# (if we ensure compatibility with AbstractAlgebra.jl)
const ExponentVector = Union{Vector}

# All supported monomial implementations in F4
const Monom = Union{ExponentVector, SparseExponentVector, PackedExponentVector}

powertype(m::Monom) = eltype(m)
powertype(::Type{M}) where {M<:Monom} = eltype(M)

#------------------------------------------------------------------------------

# Index of a monomial in the hashtable
const MonomIdx = Int32

# Hash of a monomial in the hashtable
# !!! Changing the type of hash will cause errors
const MonomHash = UInt32

# Division mask of a monomial 
const DivisionMask = UInt32

# Column index of a monomial in the matrix
const ColumnIdx = Int32
