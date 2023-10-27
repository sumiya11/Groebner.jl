using Logging

Groebner.logging_enabled() = true

@testset "logging" begin
    if Groebner.logging_enabled()
        R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
        f = [x * y + z, x * z + y]

        gb = [y^2 - z^2, x * z + y, x * y + z]
        gb1 = @test_logs Groebner.groebner(f, loglevel=-3)
        gb2 = @test_logs Groebner.groebner(f, loglevel=Int8(-3))
        prev_logger = global_logger(ConsoleLogger(stdout, Logging.Warn))
        gb3 = @test_logs Groebner.groebner(f, loglevel=Int8(-3))
        global_logger(prev_logger)
        @test gb == gb1 == gb2 == gb3
        @test_throws AssertionError Groebner.groebner(f, loglevel=:none)

        nf = Groebner.normalform(gb, f)
        nf1 = @test_logs Groebner.normalform(gb, f, loglevel=-3)
        @test nf == nf1
        @test_throws AssertionError Groebner.normalform(gb, f, loglevel=:none)

        flag1 = Groebner.isgroebner(gb)
        flag2 = @test_logs Groebner.isgroebner(gb, loglevel=-3)
        @test flag1 == flag2
        @test_throws AssertionError Groebner.isgroebner(f, loglevel=:abcd)

        f = [x^3, y^2, z]
        basis1 = Groebner.kbase(f)
        basis2 = @test_logs Groebner.kbase(f, loglevel=-3)
        @test basis1 == basis2
        @test_throws AssertionError Groebner.kbase(f, loglevel=:pkrst)
    else
        @info "Logging is disabled in Groebner.jl. Skipping tests for logging"
        @test true
    end
end
