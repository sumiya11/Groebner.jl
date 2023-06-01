
using Logging

@testset "groebner maxpairs" begin
    s = Groebner.noonn(5, ordering=:degrevlex, ground=GF(2^31 - 1))
    gb = Groebner.groebner(s)
    @test gb == Groebner.groebner(s, maxpairs=100)
    @test gb == Groebner.groebner(s, maxpairs=10)
    @test gb == Groebner.groebner(s, maxpairs=2)
    @test gb == Groebner.groebner(s, maxpairs=1)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)

    s = Groebner.katsuran(5, ordering=:degrevlex, ground=QQ)
    gb = Groebner.groebner(s)
    @test gb == Groebner.groebner(s, maxpairs=100)
    @test gb == Groebner.groebner(s, maxpairs=10)
    @test gb == Groebner.groebner(s, maxpairs=2)
    @test gb == Groebner.groebner(s, maxpairs=1)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)

    s = Groebner.katsuran(6, ordering=:deglex, ground=QQ)
    gb = Groebner.groebner(s)
    @test gb == Groebner.groebner(s, maxpairs=100)
    @test gb == Groebner.groebner(s, maxpairs=10)
    @test gb == Groebner.groebner(s, maxpairs=2)
    @test gb == Groebner.groebner(s, maxpairs=1)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)

    s = Groebner.cyclicn(5, ordering=:deglex, ground=QQ)
    gb = Groebner.groebner(s)
    @test gb == Groebner.groebner(s, maxpairs=100, linalg=:exact)
    @test gb == Groebner.groebner(s, maxpairs=10, linalg=:exact)
    @test gb == Groebner.groebner(s, maxpairs=2, linalg=:exact)
    @test gb == Groebner.groebner(s, maxpairs=2, linalg=:prob)
    @test gb == Groebner.groebner(s, maxpairs=1, linalg=:exact)
    @test gb == Groebner.groebner(s, maxpairs=1, linalg=:prob)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)

    s = Groebner.noonn(7, ordering=:degrevlex, ground=GF(2^31 - 1))
    gb = Groebner.groebner(s)
    @test gb == Groebner.groebner(s, maxpairs=100)
    @test gb == Groebner.groebner(s, maxpairs=10)
    @test gb == Groebner.groebner(s, maxpairs=1)
    @test gb == Groebner.groebner(s, maxpairs=1)

    R, (x, y, z) = QQ["x", "y", "z"]
    for i in 1:100
        s         = [x * y + z * (1 // BigInt(2)^i), x * z + y * (1 // BigInt(3)^i), y * z + x * (1 // BigInt(5)^i)]
        gblex     = Groebner.groebner(s, ordering=Groebner.Lex())
        gbdlex    = Groebner.groebner(s, ordering=Groebner.DegLex())
        gbdrevlex = Groebner.groebner(s, ordering=Groebner.DegRevLex())
        for maxpairs in [2^10, 10, 5, 1]
            @test gbdlex ==
                  Groebner.groebner(s, ordering=Groebner.DegLex(), maxpairs=maxpairs)
            @test gblex == Groebner.groebner(s, ordering=Groebner.Lex(), maxpairs=maxpairs)
            @test gbdrevlex ==
                  Groebner.groebner(s, ordering=Groebner.DegRevLex(), maxpairs=maxpairs)
        end
    end
end
