
#------------------------------------------------------------------------------

# Polynomial coefficient from finite field type.
# Among those, UInt64 is the default type used in modular computations
#=
    The most suitable type is chosen based on the size of the
    field characteristic of computation:

    2^1  <= char < 2^4   ==>   UInt8
    2^4  <= char < 2^8   ==>   UInt16
    2^8  <= char < 2^16  ==>   UInt32
    2^16 <= char < 2^32  ==>   UInt32
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

# to be done with Gleb !
# const CoeffParam = something

# All supported coefficient types in F4
const Coeff = Union{CoeffFF, CoeffZZ, CoeffQQ}

#------------------------------------------------------------------------------

# Coefficient type used in a single run of modular computation.
# Must be a subtype of CoeffFF
const CoeffModular = UInt64

@assert CoeffModular <: CoeffFF

#------------------------------------------------------------------------------

# Polynomial monomial exponent vector type;
#       works for total degrees up to 65535,
#       that can be a limitation for very large problems
#       (the problem has not been encountered yet)
const ExponentVector = Vector{UInt16}
const Degree = eltype(ExponentVector)

# Hashtable index type: is basically an index for ExponentVector in a hashtable
# get(h::Hashtable, i::ExponentIdx) --> v::ExponentVector
const ExponentIdx = Int32

const DivisionMask = UInt32

# Hash of ExponentVector
const ExponentHash = UInt32

#------------------------------------------------------------------------------

#=
    MonomsVector of a zero polynomial is an empty `MonomsVector` object.
    CoeffsVector of a zero polynomial{T} is an empty `CoeffsVector{T}` object.
=#

# Vector of polynomial monomials
const MonomsVector = Vector{ExponentIdx}
# Vector of polynomial coefficients
const CoeffsVector{T<:Coeff} = Vector{T}

#------------------------------------------------------------------------------

# Column index of a matrix
const ColumnIdx = Int32

#------------------------------------------------------------------------------