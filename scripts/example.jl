using Revise
using AbstractAlgebra, Groebner

R, (x1,x2) = polynomial_ring(QQ, ["x1","x2"], internal_ordering=:deglex);

set = [613*x1^2*x2^3 + 1413*x1*x2^2, 428*x1^3*x2 + 529*x1*x2^2]

gb = groebner(set);

isgroebner(gb)
