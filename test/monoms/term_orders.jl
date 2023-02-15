
make_ev = Groebner.make_ev
lex = (x,y) -> Groebner.monom_isless(x, y, Groebner.Lex())
dl  = (x,y) -> Groebner.monom_isless(x, y, Groebner.DegLex())
drl = (x,y) -> Groebner.monom_isless(x, y, Groebner.DegRevLex())

implementations_to_test = [
    Groebner.PowerVector{T} where {T},
    Groebner.PackedPair1{T, UInt8} where {T},
    Groebner.PackedPair2{T, UInt8} where {T},
    Groebner.PackedPair3{T, UInt8} where {T},
    # Groebner.PackedVector{T, UInt8} where {T},
]

@testset "term orders" begin
    for T in (UInt64, UInt32, UInt16)
        for EV in implementations_to_test
            a = make_ev(EV{T}, [0])
            b = make_ev(EV{T}, [0])
            @test ! lex(a, b)
            @test ! dl(a, b)
            @test ! drl(a, b)

            a = make_ev(EV{T}, [0])
            b = make_ev(EV{T}, [1])
            @test lex(a, b)
            @test dl(a, b)
            @test drl(a, b)
            @test ! lex(b, a)
            @test ! dl(b, a)
            @test ! drl(b, a)

            2 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, [1, 1])
            b = make_ev(EV{T}, [2, 0])
            @test lex(a, b)
            @test dl(a, b)
            @test drl(a, b)
            @test ! lex(b, a)
            @test ! dl(b, a)
            @test ! drl(b, a)

            2 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, [1, 1])
            b = make_ev(EV{T}, [0, 2])
            @test ! lex(a, b)
            @test ! dl(a, b)
            @test ! drl(a, b)
            @test lex(b, a)
            @test dl(b, a)
            @test drl(b, a)
            
            2 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, [1, 1])
            b = make_ev(EV{T}, [1, 1])
            @test ! lex(a, b)
            @test ! dl(a, b)
            @test ! drl(a, b)

            2 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, [1, 1])
            b = make_ev(EV{T}, [2, 2])
            @test lex(a, b)
            @test dl(a, b)
            @test drl(a, b)
            @test ! lex(b, a)
            @test ! dl(b, a)
            @test ! drl(b, a)

            3 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, [1, 0, 2])
            b = make_ev(EV{T}, [2, 0, 1])
            @test lex(a, b)
            @test dl(a, b)
            @test drl(a, b)

            5 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, ones(UInt, 5))
            b = make_ev(EV{T}, ones(UInt, 5))
            @test ! lex(a, b)
            @test ! dl(a, b)
            @test ! drl(a, b)

            25 > Groebner.capacity(EV{T}) && continue
            a = make_ev(EV{T}, ones(UInt, 25))
            b = make_ev(EV{T}, ones(UInt, 25))
            @test ! lex(a, b)
            @test ! dl(a, b)
            @test ! drl(a, b)
        end

        # test that different implementations agree
        for n in 1:2
            k = rand(1:100)
            
            implementations_to_test_local = filter(
                xxx -> Groebner.capacity(xxx{T}) >= k,
                implementations_to_test
            )

            t = div(typemax(UInt8), k) - 1
            x, y = rand(1:t, k), rand(1:t, k)
            if sum(x) >= Groebner._overflow_threshold(UInt8)
                continue
            end
            if sum(y) >= Groebner._overflow_threshold(UInt8)
                continue
            end
            as = [
                make_ev(EV{T}, x)
                for EV in implementations_to_test_local
            ]
            bs = [
                make_ev(EV{T}, y)
                for EV in implementations_to_test_local
            ]

            @test length(unique(map(lex, as, bs))) == 1
            @test length(unique(map(dl, as, bs))) == 1
            @test length(unique(map(drl, as, bs))) == 1

            # test that a < b && b < a does not happen
            for (a, b) in zip(as, bs)
                k > Groebner.capacity(a) && continue
                if lex(a, b)
                    @test ! lex(b, a)
                end
                if dl(a, b)
                    @test ! dl(b, a)
                end
                if drl(a, b)
                    @test ! drl(b, a)
                end
            end
        end
    end
end
