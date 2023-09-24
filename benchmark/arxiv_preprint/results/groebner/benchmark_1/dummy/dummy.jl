# dummy
#! format: off
using AbstractAlgebra, Groebner

ring, (x, y) = PolynomialRing(
    GF(2147483647), 
    ["x", "y"], 
    ordering=:degrevlex
)
system = [
	x^2 + y^2 + 2147483646,
	x + y
]
