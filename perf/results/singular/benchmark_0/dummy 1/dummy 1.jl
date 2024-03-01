# dummy 1
#! format: off
using AbstractAlgebra, Groebner

ring, (x, y) = polynomial_ring(
    GF(7), 
    ["x", "y"], 
    internal_ordering=:degrevlex
)
system = [
	x^2 + y^2 + 6,
	x + y
]
