
using AbstractAlgebra

@testset "ff isgroebner" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x","y","z"], ordering=:degrevlex)

    @test Groebner.isgroebner([R(1)])
    @test Groebner.isgroebner([x])
    @test Groebner.isgroebner([x, y])
    @test Groebner.isgroebner([x, x])
    @test Groebner.isgroebner([x, y, z])

    @test !Groebner.isgroebner([x + y, x])
    @test Groebner.isgroebner([z*x, z*x, R(1)])

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))

    fs = Groebner.change_ordering(Groebner.rootn(5, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test Groebner.isgroebner(Groebner.groebner(fs, reduced=false))

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))

    fs = Groebner.change_ordering(Groebner.noonn(3, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test Groebner.isgroebner(Groebner.groebner(fs, reduced=false))

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))

end
