# dummy
#! format: off
using AbstractAlgebra, Groebner

ring, (x, y) = PolynomialRing(
    QQ, 
    ["x", "y"], 
    ordering=:degrevlex
)
system = [
	x^2 + y^2 - 1,
	x + y
]
