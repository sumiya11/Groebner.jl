# How do we test the output printed to stdout?..

function test_timings()
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    f = [x * y + z, x * z + y]

    gb = [y^2 - z^2, x * z + y, x * y + z]
    gb1 = Groebner.groebner(f, statistics=:no)
    gb2 = Groebner.groebner(f, statistics=:timings)
    gb3 = Groebner.groebner(f, statistics=:all)
    @test gb == gb1 == gb2 == gb3
    @test_throws AssertionError Groebner.groebner(f, statistics=:alll)

    nf = Groebner.normalform(gb, f)
    nf1 = Groebner.normalform(gb, f, statistics=:no)
    nf2 = Groebner.normalform(gb, f, statistics=:timings)
    nf3 = Groebner.normalform(gb, f, statistics=:all)
    @test nf == nf1 == nf2 == nf3
    @test_throws AssertionError Groebner.normalform(gb, f, statistics=:alll)

    flag1 = Groebner.isgroebner(gb, statistics=:no)
    flag2 = Groebner.isgroebner(gb, statistics=:timings)
    flag3 = Groebner.isgroebner(gb, statistics=:all)
    @test flag1 == flag2 == flag3
    @test_throws AssertionError Groebner.isgroebner(f, statistics=:alll)

    f = [x^3, y^2, z]
    basis1 = Groebner.kbase(f, statistics=:no)
    basis2 = Groebner.kbase(f, statistics=:timings)
    basis3 = Groebner.kbase(f, statistics=:all)
    @test basis1 == basis2 == basis3
    @test_throws AssertionError Groebner.kbase(f, statistics=:alll)
end

@testset "collect statistics: timings" begin
    Groebner.performance_counters_enabled() = true
    test_timings()
    Groebner.performance_counters_enabled() = false
    test_timings()
end
