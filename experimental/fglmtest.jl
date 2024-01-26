using AbstractAlgebra, Groebner

R, (x, y, z, t) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x", "y", "z", "t"])
sys = [
    y^2 * z + 2 * x * y * t - 2 * x - z,
    -x^3 * z + 4 * x * y^2 * z + 4 * x^2 * y * t + 2 * y^3 * t + 4 * x^2 - 10 * y^2 +
    4 * x * z - 10 * y * t + 2,
    2 * y * z * t + x * t^2 - x - 2 * z,
    -x * z^3 + 4 * y * z^2 * t + 4 * x * z * t^2 + 2 * y * t^3 + 4 * x * z + 4 * z^2 -
    10 * y * t - 10 * t^2 + 2
]

gb_lex = Groebner.groebner(sys, ordering=Groebner.Lex())
gb_drl = Groebner.groebner(sys, ordering=Groebner.DegRevLex())
gb_fglm = Groebner.fglm(gb_drl, Groebner.DegRevLex(), Groebner.Lex())

@info "" gb_fglm
