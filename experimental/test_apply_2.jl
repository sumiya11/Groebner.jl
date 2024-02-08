using AbstractAlgebra
include("../src/Groebner.jl")

R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)
p, p2 = 3, 19

sys = [-11 * x * y + 53 * y, 83 * x * y + x - 70 * y]
sys_mod_p = map(
    f -> map_coefficients(c -> GF(p)(numerator(c)), f),
    [x * y + 2 * y, 2 * x * y + x + 2 * y]
)
sys_mod_p2 = map(
    f -> map_coefficients(c -> GF(p2)(numerator(c)), f),
    [8 * x * y + 15 * y, 7 * x * y + x + 6 * y]
)

trace, gb = Groebner.groebner_learn(sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p)
@enter flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2)
