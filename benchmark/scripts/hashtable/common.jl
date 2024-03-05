
const Monom = Vector{UInt8}

# The idenfifier of a monomial. This idenfifier is guaranteed to be unique
# within a particular hashtable. This allows one to use this idenfifier when
# working with monomials
const MonomId = Int32
