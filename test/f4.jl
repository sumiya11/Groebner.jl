
using AbstractAlgebra

#------------------------------------------------------------------------------

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
#=
    ans is
    y^2 + x,
    x^2 + x,
    x*y + x
=#


fs = [
    x + y,
    x^2 + y
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x + y,
    x^2 + y
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y,
    x*y + x
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    5y,
    4x*y + 8x
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@assert gb == [R(1)]

fs = [
    x*y + y,
    y - 5
]
gb = GroebnerBases.f4(fs)
println(gb)
@assert GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@assert gb == [R(1)]

#------------------------------------------------------------------------------
