# cyclic 8
#! format: off
using AbstractAlgebra, Groebner

ring, (z1, z2, z3, z4, z5, z6, z7, z8) = PolynomialRing(
    GF(2147483647), 
    ["z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8"], 
    ordering=:degrevlex
)
system = [
	z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8,
	z1*z2 + z1*z8 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z8,
	z1*z2*z3 + z1*z2*z8 + z1*z7*z8 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7 + z6*z7*z8,
	z1*z2*z3*z4 + z1*z2*z3*z8 + z1*z2*z7*z8 + z1*z6*z7*z8 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7 + z5*z6*z7*z8,
	z1*z2*z3*z4*z5 + z1*z2*z3*z4*z8 + z1*z2*z3*z7*z8 + z1*z2*z6*z7*z8 + z1*z5*z6*z7*z8 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7 + z4*z5*z6*z7*z8,
	z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z8 + z1*z2*z3*z4*z7*z8 + z1*z2*z3*z6*z7*z8 + z1*z2*z5*z6*z7*z8 + z1*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z8,
	z1*z2*z3*z4*z5*z6*z7 + z1*z2*z3*z4*z5*z6*z8 + z1*z2*z3*z4*z5*z7*z8 + z1*z2*z3*z4*z6*z7*z8 + z1*z2*z3*z5*z6*z7*z8 + z1*z2*z4*z5*z6*z7*z8 + z1*z3*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7*z8,
	z1*z2*z3*z4*z5*z6*z7*z8 + 2147483646
]
