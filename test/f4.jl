
using AbstractAlgebra
using Logging

global_logger(ConsoleLogger(stderr, Logging.Error))

@testset "f4" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
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
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x + y,
    x^2 + y
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y,
    x*y + x
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    5y,
    4x*y + 8x
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@test gb == [R(1)]

fs = [
    x*y + y,
    y - 5
]
gb = GroebnerBases.groebner(fs)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

root = GroebnerBases.change_ordering(GroebnerBases.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

root = GroebnerBases.change_ordering(GroebnerBases.rootn(4, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

for proot in Combinatorics.permutations(root)
    gb = GroebnerBases.groebner(proot)
    @test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
end

root = GroebnerBases.change_ordering(GroebnerBases.rootn(8, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

noon = GroebnerBases.change_ordering(GroebnerBases.noon3(ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(noon)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

for pnoon in Combinatorics.permutations(noon)
    gb = GroebnerBases.groebner(pnoon)
    @test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
end

noon = GroebnerBases.change_ordering(GroebnerBases.noon4(ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(noon)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

end
