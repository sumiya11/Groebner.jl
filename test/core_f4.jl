
using AbstractAlgebra


@testset "ff f4 noreduce degrevlex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
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
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

fs = [
    x + y,
    x^2 + y
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

fs = [
    y,
    x*y + x
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

fs = [
    5y,
    4x*y + 8x
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
@test gb == [R(1)]

fs = [
    x*y + y,
    y - 5
]
gb = Groebner.groebner(fs, reduced=false)
# println(gb)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

root = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(root, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

root = Groebner.change_ordering(Groebner.rootn(4, ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(root, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

for proot in Groebner.Combinatorics.permutations(root)
    gb = Groebner.groebner(proot, reduced=false)
    @test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
end

root = Groebner.change_ordering(Groebner.rootn(8, ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(root, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

noon = Groebner.change_ordering(Groebner.noonn(3, ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(noon, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

for pnoon in Groebner.Combinatorics.permutations(noon)
    gb = Groebner.groebner(pnoon, reduced=false)
    @test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
end

noon = Groebner.change_ordering(Groebner.noonn(4, ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(noon, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

eco = Groebner.change_ordering(Groebner.eco5(), :degrevlex)
gb = Groebner.groebner(eco, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

ku = Groebner.change_ordering(Groebner.ku10(ground=GF(2^31 - 1)), :degrevlex)
gb = Groebner.groebner(ku, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

end


@testset "ff f4 noreduce lex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = Groebner.groebner(fs, reduced=false)
@test gb ≂ [x + y^2, x*y + 2147483646*y^2, y^3 + y^2]

fs = [
    x + y,
    x^2 + y
]
gb = Groebner.groebner(fs, reduced=false)
@test gb ≂ [x + y, x^2 + y, y^2 + y]

fs = [
    y,
    x*y + x
]
gb = Groebner.groebner(fs, reduced=false)
@test gb ≂ [x, y]

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = Groebner.groebner(fs, reduced=false)
@test gb ≂ [ y^2 + 1073741825, x^2 + 5]

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
@test gb == [R(1)]

fs = [
    x^2*y^2,
    2*x*y^2 + 3*x*y
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
@test gb ≂ [ x*y^2 + 1073741825*x*y, x^2*y]


root = Groebner.rootn(3, ground=GF(2^31 - 1))
gb = Groebner.groebner(root, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

root = Groebner.rootn(4, ground=GF(2^31 - 1))
gb = Groebner.groebner(root, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

#=
noon = Groebner.noon3(ground=GF(2^31 - 1))
gb = Groebner.groebner(noon, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
=#

end

@testset "ff f4 noreduce deglex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:deglex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = Groebner.groebner(fs, reduced=true)
@test Groebner._isgroebner_reference(gb, initial_gens=fs)

fs = [
    x + y,
    x^2 + y
]
gb = Groebner.groebner(fs, reduced=true)
@test Groebner._isgroebner_reference(gb, initial_gens=fs)

fs = [
    y,
    x*y + x
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(gb, initial_gens=fs)

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(gb, initial_gens=fs)

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
@test Groebner._isgroebner_reference(gb, initial_gens=fs)

fs = [
    x^2*y^2,
    2*x*y^2 + 3*x*y
]
gb = Groebner.groebner(fs, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))
@test Groebner._isgroebner_reference(gb, initial_gens=fs)


root = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31 - 1)), :deglex)
gb = Groebner.groebner(root, reduced=true)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))


root = Groebner.change_ordering(Groebner.rootn(4, ground=GF(2^31 - 1)), :deglex)
gb = Groebner.groebner(root, reduced=true)
@test Groebner._isgroebner_reference(gb, initial_gens=root)

noon = Groebner.change_ordering(Groebner.noonn(3, ground=GF(2^31 - 1)), :deglex)
gb = Groebner.groebner(noon, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

noon = Groebner.change_ordering(Groebner.noonn(4, ground=GF(2^31 - 1)), :deglex)
gb = Groebner.groebner(noon, reduced=false)
@test Groebner._isgroebner_reference(Groebner._reducegb_reference(gb))

end

@testset "ff f4 corner cases" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

    fs = [x, x, x]
    @test Groebner.groebner(fs) ≂ [x]

    fs = [x, x, y, y, y, x, x, y]
    @test Groebner.groebner(fs) ≂ [x, y]

    fs = [R(1)]
    @test Groebner.groebner(fs) ≂ [R(1)]

    fs = [x^2 + y, y*x - 1, R(1), y^4]
    @test Groebner.groebner(fs) ≂ [R(1)]

    fs = [2x, 3y, 4x + 5y]
    @test Groebner.groebner(fs) ≂ [x, y]


    R, (x1, x2) = PolynomialRing(GF(2^31 - 1), ["x1", "x2"], ordering=:degrevlex)

    fs = [1*x1^2*x2^2 + 2*x1^2*x2,
        1*x1^2*x2^2 + 3*x1^2*x2 + 5*x1*x2^2,
        1*x1^2*x2^2]
    gb = Groebner.groebner(fs)
    @test gb ≂ [x1*x2^2, x1^2*x2]

    fs = [x1*x2^2 + x1*x2, x1^2*x2^2 + 2*x1*x2, 2*x1^2*x2^2 + x1*x2]
    gb = Groebner.groebner(fs)
    @test gb ≂ [x1*x2]

end
