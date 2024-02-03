# jason210
#! format: off
using AbstractAlgebra, Groebner

ring, (x1, x2, x3, x4, x5, x6, x7, x8) = polynomial_ring(
    GF(1073741827), 
    ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"], 
    ordering=:degrevlex
)
system = [
	x1^2*x3^4 + x2^2*x4^4 + x1*x2*x3^2*x5^2 + x1*x2*x4^2*x6^2 + x1*x2*x3*x4*x5*x7 + x1*x2*x3*x4*x6*x8,
	x2^6,
	x1^6
]

