
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
end
