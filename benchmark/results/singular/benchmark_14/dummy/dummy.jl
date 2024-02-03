# dummy
#! format: off
using AbstractAlgebra, Groebner

ring, (x, y) = polynomial_ring(
    GF(1073741827), 
    ["x", "y"], 
    ordering=:degrevlex
)
system = [
	x^2 + y^2 + 1073741826,
	x + y
]

