

# Polynomial coefficient from finite field type
const CoeffFF = UInt64

# Polynomial coefficient from integer ring type
const CoeffZZ = BigInt

# Polynomial coefficient from rational field type
const CoeffQQ = Rational{BigInt}

# All supported coefficient types in F4
const Coeff = Union{CoeffFF, CoeffZZ, CoeffQQ}

# Polynomial term exponent vector type
#       works for total degrees up to 65536,
#       that can be a limitation for very large problems
const ExponentVector = Vector{UInt16}
