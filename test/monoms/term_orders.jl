
make_ev = Groebner.make_ev
lex = (x,y) -> Groebner.monom_isless(x, y, Groebner.Lex())
dl  = (x,y) -> Groebner.monom_isless(x, y, Groebner.DegLex())
drl = (x,y) -> Groebner.monom_isless(x, y, Groebner.DegRevLex())

implementations_to_test = [
    Groebner.PowerVector{T} where {T},
    Groebner.PackedPair1{T, UInt8} where {T},
    Groebner.PackedPair2{T, UInt8} where {T},
    Groebner.PackedPair3{T, UInt8} where {T},
]

@testset "term orders: Lex, DegLex, DegRevLex" begin
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

function test_orderings(v1, v2, ords_to_test)
    for wo in ords_to_test
        ord = wo.ord
        ans = wo.ans
        @test Groebner.monom_isless(v1, v2, ord) == ans[1]
        @test Groebner.monom_isless(v2, v1, ord) == ans[2]
    end
end

@testset "term orders: WeightedOrdering" begin    
    pv = Groebner.PowerVector{T} where {T}

    v1 = Groebner.make_ev(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.make_ev(pv{UInt64}, [3, 2, 1])

    ords_to_test = [
        (ord=Groebner.WeightedOrdering([1, 1, 1]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([0, 0, 1]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([0, 1, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 1, 5]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
    ]

    test_orderings(v1, v2, ords_to_test)

    v1 = Groebner.make_ev(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.make_ev(pv{UInt64}, [1, 2, 3])

    ords_to_test = [
        (ord=Groebner.WeightedOrdering([1, 1, 1]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([0, 0, 1]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([0, 1, 0]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 1, 5]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[false, false]),
    ]

    test_orderings(v1, v2, ords_to_test)

    v1 = Groebner.make_ev(pv{UInt64}, [1, 2, 3, 0, 0, 7])
    v2 = Groebner.make_ev(pv{UInt64}, [4, 5, 0, 0, 1, 4])

    ords_to_test = [
        (ord=Groebner.WeightedOrdering([1, 1, 1]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([0, 0, 1]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([0, 1, 0]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 1, 5]),ans=[false, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[false, false]),
    ]

    test_orderings(v1, v2, ords_to_test)

end

@testset "term orders: BlockOrdering" begin
    pv = Groebner.PowerVector{T} where {T}

    v1 = Groebner.make_ev(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.make_ev(pv{UInt64}, [3, 2, 1])

    wo_to_test = [
        (ord=Groebner.WeightedOrdering([1, 1, 1]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([0, 0, 1]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([0, 1, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 1, 5]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
    ]
    for wo in wo_to_test
        ord = wo.ord
        ans = wo.ans
        @test Groebner.monom_isless(v1, v2, ord) == ans[1]
        @test Groebner.monom_isless(v2, v1, ord) == ans[2]
    end

end


@testset "term orders: MatrixOrdering" begin
    pv = Groebner.PowerVector{T} where {T}

    v1 = Groebner.make_ev(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.make_ev(pv{UInt64}, [3, 2, 1])

    wo_to_test = [
        (ord=Groebner.WeightedOrdering([1, 1, 1]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([0, 0, 1]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([0, 1, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
        (ord=Groebner.WeightedOrdering([1, 1, 5]),ans=[false, true]),
        (ord=Groebner.WeightedOrdering([1, 0, 0]),ans=[true, false]),
    ]
    for wo in wo_to_test
        ord = wo.ord
        ans = wo.ans
        @test Groebner.monom_isless(v1, v2, ord) == ans[1]
        @test Groebner.monom_isless(v2, v1, ord) == ans[2]
    end

end
