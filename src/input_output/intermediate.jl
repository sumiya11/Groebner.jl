# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Intermediate polynomial representation (ir)

# Polynomials in the intermediate representation are represented 
# by coefficients (UInt64) and exponent vectors (Vector{UInt64}).

"""
    PolyRing

A polynomial ring.

## Example

4 variables, modulo 2^30 + 3, degrevlex order.

```julia
ring = Groebner.PolyRing(4, Groebner.DegRevLex(), 2^30+3)
```

4 variables, the rationals, a block order.

```julia
ord = Groebner.DegRevLex([1,2]) * Groebner.DegRevLex([3,4])
ring = Groebner.PolyRing(4, ord, 0)
```
"""
mutable struct PolyRing{
    Ord <: AbstractMonomialOrdering,
    C <: Union{CoeffZp, CompositeCoeffZp}
}
    nvars::Int
    ord::Ord
    ch::C
end
