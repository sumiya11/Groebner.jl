# This file was generated, do not modify it. # hide
using AbstractAlgebra
_, (t, x, y) = PolynomialRing(QQ, ["t", "x", "y"], ordering=:lex)

groebner([t^2*y - 2t + y, t^2*x + t^2 + x - 1])