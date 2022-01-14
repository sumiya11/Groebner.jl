
using AbstractAlgebra


@testset "ff f4 noreduce degrevlex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
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
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

fs = [
    x + y,
    x^2 + y
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

fs = [
    y,
    x*y + x
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

fs = [
    5y,
    4x*y + 8x
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
@test gb == [R(1)]

fs = [
    x*y + y,
    y - 5
]
gb = FastGroebner.groebner(fs, reduced=false)
# println(gb)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

root = FastGroebner.change_ordering(FastGroebner.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(root, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

root = FastGroebner.change_ordering(FastGroebner.rootn(4, ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(root, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

for proot in FastGroebner.Combinatorics.permutations(root)
    gb = FastGroebner.groebner(proot, reduced=false)
    @test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
end

root = FastGroebner.change_ordering(FastGroebner.rootn(8, ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(root, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

noon = FastGroebner.change_ordering(FastGroebner.noon3(ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(noon, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

for pnoon in FastGroebner.Combinatorics.permutations(noon)
    gb = FastGroebner.groebner(pnoon, reduced=false)
    @test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
end

noon = FastGroebner.change_ordering(FastGroebner.noon4(ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(noon, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

eco = FastGroebner.change_ordering(FastGroebner.eco5(), :degrevlex)
gb = FastGroebner.groebner(eco, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

ku = FastGroebner.change_ordering(FastGroebner.ku10(ground=GF(2^31 - 1)), :degrevlex)
gb = FastGroebner.groebner(ku, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

end


@testset "ff f4 noreduce lex" begin

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

fs = [
    x + y^2,
    x*y - y^2
]
gb = FastGroebner.groebner(fs, reduced=false)
@test gb ≂ [x + y^2, x*y + 2147483646*y^2, y^3 + y^2]

fs = [
    x + y,
    x^2 + y
]
gb = FastGroebner.groebner(fs, reduced=false)
@test gb ≂ [x + y, x^2 + y, y^2 + y]

fs = [
    y,
    x*y + x
]
gb = FastGroebner.groebner(fs, reduced=false)
@test gb ≂ [x, y]

fs = [
    x^2 + 5,
    2y^2 + 3
]
gb = FastGroebner.groebner(fs, reduced=false)
@test gb ≂ [ y^2 + 1073741825, x^2 + 5]

fs = [
    y^2 + x,
    x^2*y + y,
    y + 1
]
gb = FastGroebner.groebner(fs, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
@test gb == [R(1)]

fs = [
    x^2*y^2,
    2*x*y^2 + 3*x*y
]
gb = FastGroebner.groebner(fs, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
@test gb ≂ [ x*y^2 + 1073741825*x*y, x^2*y]


root = FastGroebner.rootn(3, ground=GF(2^31 - 1))
gb = FastGroebner.groebner(root, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

root = FastGroebner.rootn(4, ground=GF(2^31 - 1))
gb = FastGroebner.groebner(root, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))

#=
noon = FastGroebner.noon3(ground=GF(2^31 - 1))
gb = FastGroebner.groebner(noon, reduced=false)
@test FastGroebner.isgroebner(FastGroebner.reducegb(gb))
=#

end

@testset "ff f4 corner cases" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

    fs = [x, x, x]
    @test FastGroebner.groebner(fs) ≂ [x]

    fs = [x, x, y, y, y, x, x, y]
    @test FastGroebner.groebner(fs) ≂ [x, y]

    fs = [R(1)]
    @test FastGroebner.groebner(fs) ≂ [R(1)]

    fs = [x^2 + y, y*x - 1, R(1), y^4]
    @test FastGroebner.groebner(fs) ≂ [R(1)]

    fs = [2x, 3y, 4x + 5y]
    @test FastGroebner.groebner(fs) ≂ [x, y]


    R, (x1, x2) = PolynomialRing(GF(2^31 - 1), ["x1", "x2"], ordering=:degrevlex)

    fs = [1*x1^2*x2^2 + 2*x1^2*x2,
        1*x1^2*x2^2 + 3*x1^2*x2 + 5*x1*x2^2,
        1*x1^2*x2^2]
    gb = FastGroebner.groebner(fs)
    @test gb ≂ [x1*x2^2, x1^2*x2]

    fs = [x1*x2^2 + x1*x2, x1^2*x2^2 + 2*x1*x2, 2*x1^2*x2^2 + x1*x2]
    gb = FastGroebner.groebner(fs)
    @test gb ≂ [x1*x2]

end
