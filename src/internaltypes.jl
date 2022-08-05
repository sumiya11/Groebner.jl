
#------------------------------------------------------------------------------

# Polynomial coefficient from finite field type.
# Among those, UInt64 is the default type used in modular computations
const CoeffFF = Union{UInt8, UInt16, UInt32, UInt64}

# Polynomial coefficient from integer ring type
const CoeffZZ = Union{BigInt}

# Polynomial coefficient from rational field type
const CoeffQQ = Union{Rational{BigInt}}

# to be done with Gleb !
# const CoeffParam = something

# All supported coefficient types in F4
const Coeff = Union{CoeffFF, CoeffZZ, CoeffQQ}

# Field characteristic
const FieldCharacteristic = UInt64

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

# Vector of polynomial monomials
const MonomsVector = Vector{ExponentIdx}
# Vector of polynomial coefficients
const CoeffsVector{T<:Coeff} = Vector{T}

#------------------------------------------------------------------------------

# Column index of a matrix
const ColumnIdx = Int32

#------------------------------------------------------------------------------