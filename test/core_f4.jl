
using AbstractAlgebra


@testset "ff f4 noreduce degrevlex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = GroebnerBases.groebner(fs, reduced=false)
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
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x + y,
    x^2 + y
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y,
    x*y + x
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    5y,
    4x*y + 8x
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@test gb == [R(1)]

fs = [
    x*y + y,
    y - 5
]
gb = GroebnerBases.groebner(fs, reduced=false)
# println(gb)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

root = GroebnerBases.change_ordering(GroebnerBases.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

root = GroebnerBases.change_ordering(GroebnerBases.rootn(4, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

for proot in GroebnerBases.Combinatorics.permutations(root)
    gb = GroebnerBases.groebner(proot, reduced=false)
    @test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
end

root = GroebnerBases.change_ordering(GroebnerBases.rootn(8, ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(root, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

noon = GroebnerBases.change_ordering(GroebnerBases.noon3(ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(noon, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

for pnoon in GroebnerBases.Combinatorics.permutations(noon)
    gb = GroebnerBases.groebner(pnoon, reduced=false)
    @test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
end

noon = GroebnerBases.change_ordering(GroebnerBases.noon4(ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(noon, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

eco = GroebnerBases.change_ordering(GroebnerBases.eco5(), :degrevlex)
gb = GroebnerBases.groebner(eco, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

ku = GroebnerBases.change_ordering(GroebnerBases.ku10(ground=GF(2^31 - 1)), :degrevlex)
gb = GroebnerBases.groebner(ku, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

end


@testset "ff f4 noreduce lex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test gb ≂ [x + y^2, x*y + 2147483646*y^2, y^3 + y^2]

fs = [
    x + y,
    x^2 + y
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test gb ≂ [x + y, x^2 + y, y^2 + y]

fs = [
    y,
    x*y + x
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test gb ≂ [x, y]

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test gb ≂ [ y^2 + 1073741825, x^2 + 5]

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@test gb == [R(1)]

fs = [
    x^2*y^2,
    2*x*y^2 + 3*x*y
]
gb = GroebnerBases.groebner(fs, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
@test gb ≂ [ x*y^2 + 1073741825*x*y, x^2*y]


root = GroebnerBases.rootn(3, ground=GF(2^31 - 1))
gb = GroebnerBases.groebner(root, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

root = GroebnerBases.rootn(4, ground=GF(2^31 - 1))
gb = GroebnerBases.groebner(root, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))

#=
noon = GroebnerBases.noon3(ground=GF(2^31 - 1))
gb = GroebnerBases.groebner(noon, reduced=false)
@test GroebnerBases.isgroebner(GroebnerBases.reducegb(gb))
=#

end

@testset "ff f4 corner cases" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

    fs = [x, x, x]
    @test GroebnerBases.groebner(fs) ≂ [x]

    fs = [x, x, y, y, y, x, x, y]
    @test GroebnerBases.groebner(fs) ≂ [x, y]

    fs = [R(1)]
    @test GroebnerBases.groebner(fs) ≂ [R(1)]

    fs = [x^2 + y, y*x - 1, R(1), y^4]
    @test GroebnerBases.groebner(fs) ≂ [R(1)]

    fs = [2x, 3y, 4x + 5y]
    @test GroebnerBases.groebner(fs) ≂ [x, y]


    R, (x1, x2) = PolynomialRing(GF(2^31 - 1), ["x1", "x2"], ordering=:degrevlex)

    fs = [1*x1^2*x2^2 + 2*x1^2*x2,
        1*x1^2*x2^2 + 3*x1^2*x2 + 5*x1*x2^2,
        1*x1^2*x2^2]
    gb = GroebnerBases.groebner(fs)
    @test gb ≂ [x1*x2^2, x1^2*x2]

    # TODO: reduction does not work
    fs = [x1*x2^2 + x1*x2, x1^2*x2^2 + 2*x1*x2, 2*x1^2*x2^2 + x1*x2]
    gb = GroebnerBases.groebner(fs)
    @test gb ≂ [x1*x2]


end
