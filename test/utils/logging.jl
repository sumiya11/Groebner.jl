using Logging

@testset "logging" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    f = [x * y + z, x * z + y]
    gb = [y^2 - z^2, x * z + y, x * y + z]

    # Simple tests with logging enabled
    Groebner.logging_enabled() = true
    if Groebner.logging_enabled()
        gb1 = @test_logs Groebner.groebner(f, loglevel=-3)
        prev_logger = global_logger(ConsoleLogger(stdout, Logging.Warn))
        gb2 = @test_logs Groebner.groebner(f, loglevel=Int8(-3))
        global_logger(prev_logger)
        @test gb == gb1 == gb2
        @test_throws AssertionError Groebner.groebner(f, loglevel=:none)

        nf = Groebner.normalform(gb, f)
        nf1 = @test_logs Groebner.normalform(gb, f, loglevel=-3)
        @test nf == nf1
        @test_throws AssertionError Groebner.normalform(gb, f, loglevel=:none)

        flag1 = Groebner.isgroebner(gb)
        flag2 = @test_logs Groebner.isgroebner(gb, loglevel=-3)
        @test flag1 == flag2
        @test_throws AssertionError Groebner.isgroebner(f, loglevel=:abcd)

        gb1 = @test_logs Groebner.groebner(f, loglevel=:debug)
        gb1 = @test_logs Groebner.groebner(f, loglevel=:info)
        gb1 = @test_logs Groebner.groebner(f, loglevel=:warn)
    end

    # Simple tests with logging disabled
    Groebner.logging_enabled() = false
    if !Groebner.logging_enabled()
        gb1 = Groebner.groebner(f, loglevel=-3)
        prev_logger = global_logger(ConsoleLogger(stdout, Logging.Warn))
        gb2 = @test_nowarn Groebner.groebner(f, loglevel=Int8(-3))
        global_logger(prev_logger)
        @test gb == gb1 == gb2
        @test_throws AssertionError Groebner.groebner(f, loglevel=:none)

        nf = Groebner.normalform(gb, f)
        nf1 = Groebner.normalform(gb, f, loglevel=-3)
        @test nf == nf1
        @test_throws AssertionError Groebner.normalform(gb, f, loglevel=:none)

        flag1 = Groebner.isgroebner(gb)
        flag2 = Groebner.isgroebner(gb, loglevel=-3)
        @test flag1 == flag2
        @test_throws AssertionError Groebner.isgroebner(f, loglevel=:abcd)
    end
end
