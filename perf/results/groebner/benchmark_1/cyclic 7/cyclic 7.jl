# cyclic 7
#! format: off
using AbstractAlgebra, Groebner

ring, (z1, z2, z3, z4, z5, z6, z7) = PolynomialRing(
    GF(1073741827), 
    ["z1", "z2", "z3", "z4", "z5", "z6", "z7"], 
    internal_ordering=:degrevlex
)
system = [
	z1 + z2 + z3 + z4 + z5 + z6 + z7,
	z1*z2 + z1*z7 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7,
	z1*z2*z3 + z1*z2*z7 + z1*z6*z7 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7,
	z1*z2*z3*z4 + z1*z2*z3*z7 + z1*z2*z6*z7 + z1*z5*z6*z7 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7,
	z1*z2*z3*z4*z5 + z1*z2*z3*z4*z7 + z1*z2*z3*z6*z7 + z1*z2*z5*z6*z7 + z1*z4*z5*z6*z7 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7,
	z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z7 + z1*z2*z3*z4*z6*z7 + z1*z2*z3*z5*z6*z7 + z1*z2*z4*z5*z6*z7 + z1*z3*z4*z5*z6*z7 + z2*z3*z4*z5*z6*z7,
	z1*z2*z3*z4*z5*z6*z7 + 1073741826
]
