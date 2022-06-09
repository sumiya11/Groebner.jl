
using AbstractAlgebra

@testset "isgroebner certify=true" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x","y","z"], ordering=:degrevlex)

    @test Groebner.isgroebner([R(1)], certify=true)
    @test Groebner.isgroebner([x], certify=true)
    @test Groebner.isgroebner([x, y], certify=true)
    @test Groebner.isgroebner([x, x], certify=true)
    @test Groebner.isgroebner([x, y, z], certify=true)

    @test !Groebner.isgroebner([x + y, x], certify=true)
    @test Groebner.isgroebner([z*x, z*x, R(1)], certify=true)

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=true)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=true)
end

@testset "isgroebner certify=false" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x","y","z"], ordering=:degrevlex)

    @test Groebner.isgroebner([R(1)], certify=false)
    @test Groebner.isgroebner([x], certify=false)
    @test Groebner.isgroebner([x, y], certify=false)
    @test Groebner.isgroebner([x, x], certify=false)
    @test Groebner.isgroebner([x, y, z], certify=false)

    @test !Groebner.isgroebner([x + y, x], certify=false)
    @test Groebner.isgroebner([z*x, z*x, R(1)], certify=false)

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31-1)), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=QQ), :degrevlex)
    @test !Groebner.isgroebner(fs, certify=false)
    @test Groebner.isgroebner(Groebner.groebner(fs), certify=false)
end
