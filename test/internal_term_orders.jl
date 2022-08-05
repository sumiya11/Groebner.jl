
EV = Groebner.ExponentVector
lex = Groebner.exponent_isless_lex
dl  = Groebner.exponent_isless_dl
drl = Groebner.exponent_isless_drl

@testset "Correctness of lex, deglex, degrevlex" begin

    a = EV([0, 0])
    b = EV([0, 0])
    @test ! lex(a, b)
    @test ! dl(a, b)
    @test ! drl(a, b)

    a = EV([0, 0])
    b = EV([1, 1])
    @test lex(a, b)
    @test dl(a, b)
    @test drl(a, b)
    @test ! lex(b, a)
    @test ! dl(b, a)
    @test ! drl(b, a)

    a = EV([1, 1, 2])
    b = EV([2, 0, 2])
    @test lex(a, b)
    @test dl(a, b)
    @test drl(a, b)
    @test ! lex(b, a)
    @test ! dl(b, a)
    @test ! drl(b, a)

    a = EV([1, 1, 2])
    b = EV([0, 2, 2])
    @test ! lex(a, b)
    @test ! dl(a, b)
    @test ! drl(a, b)
    @test lex(b, a)
    @test dl(b, a)
    @test drl(b, a)
    
    a = EV([1, 1, 2])
    b = EV([1, 1, 2])
    @test ! lex(a, b)
    @test ! dl(a, b)
    @test ! drl(a, b)

    a = EV([1, 1, 2])
    b = EV([2, 2, 2])
    @test lex(a, b)
    @test dl(a, b)
    @test ! drl(a, b)
    @test ! lex(b, a)
    @test ! dl(b, a)
    @test drl(b, a)
end