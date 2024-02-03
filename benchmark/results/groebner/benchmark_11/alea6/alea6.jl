# alea6
#! format: off
using AbstractAlgebra, Groebner

ring, (x0, x1, x2, x3, x4, x5) = polynomial_ring(
    QQ, 
    ["x0", "x1", "x2", "x3", "x4", "x5"], 
    ordering=:degrevlex
)
system = [
	5*x0^2*x3 + 37*x1*x3*x4 + 32*x1*x3*x5 + 21*x3*x5 + 55*x4*x5,
	23*x1^2*x4 + 57*x1*x2*x4 + 56*x1*x4^2 + 39*x0*x1*x5 + 52*x3*x4*x5 + 10*x2^2,
	33*x0^2*x3 + 32*x1*x3^2 + 51*x1^2*x4 + 42*x0*x3*x5 + x5^3 + 51*x0^2,
	44*x0*x3^2 + 47*x1*x4^2 + 43*x3*x4^2 + 2*x2*x4*x5 + 42*x1*x3 + 12*x2*x3,
	49*x0^2*x2 + 11*x0*x1*x2 + 45*x1^2*x4 + 83*x0*x3*x4 + 54*x0*x3,
	48*x0*x2*x3 + 2*x2^2*x3 + 36*x3^3 + 59*x2^2*x5 + 17*x2 + 45*x4
]

