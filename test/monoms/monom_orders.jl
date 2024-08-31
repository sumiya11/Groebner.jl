using AbstractAlgebra

monom_construct_from_vector = Groebner.monom_construct_from_vector
lex = (x, y) -> Groebner.monom_isless(x, y, Groebner.Lex())
dl = (x, y) -> Groebner.monom_isless(x, y, Groebner.DegLex())
drl = (x, y) -> Groebner.monom_isless(x, y, Groebner.DegRevLex())

implementations_to_test = [
    Groebner.ExponentVector{T} where {T},
    Groebner.PackedTuple1{T, UInt8} where {T},
    Groebner.PackedTuple2{T, UInt8} where {T},
    Groebner.PackedTuple3{T, UInt8} where {T},
    Groebner.PackedTuple4{T, UInt8} where {T}
]

@testset "monom orders: Lex, DegLex, DegRevLex" begin
    x = [1, 1, 1, 1, 1, 0, 1, 1, 0, 1]
    y = [1, 1, 1, 1, 1, 0, 1, 1, 0, 1]
    a = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, x)
    b = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, y)
    @test !lex(a, b) && !dl(a, b) && !drl(a, b)

    x = [1, 1, 1, 1, 1, 0, 1, 2, 0, 1]
    y = [1, 1, 1, 1, 1, 0, 2, 1, 0, 1]
    a = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, x)
    b = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, y)
    @test lex(a, b) && dl(a, b) && drl(a, b)

    n = 12
    for i in 2:n
        xx = ones(Int32, n)
        xx[i] = xx[i] + 1
        yy = ones(Int32, n)
        yy[i - 1] = yy[i - 1] + 1
        a = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, xx)
        b = monom_construct_from_vector(Groebner.ExponentVector{UInt32}, yy)
        @test lex(a, b) && dl(a, b) && drl(a, b)
    end

    for T in (UInt64, UInt32, UInt16, UInt8)
        for EV in implementations_to_test
            if EV{T} <: Groebner.AbstractPackedTuple{T, UInt8} && T == UInt8
                continue
            end

            a = monom_construct_from_vector(EV{T}, [0])
            b = monom_construct_from_vector(EV{T}, [0])
            @test !drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test !lex(a, b)
                @test !dl(a, b)
            end

            a = monom_construct_from_vector(EV{T}, [0])
            b = monom_construct_from_vector(EV{T}, [1])
            @test drl(a, b)
            @test !drl(b, a)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test lex(a, b)
                @test dl(a, b)
                @test !lex(b, a)
                @test !dl(b, a)
            end

            2 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 1])
            b = monom_construct_from_vector(EV{T}, [2, 0])
            @test drl(a, b)
            @test !drl(b, a)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test lex(a, b)
                @test dl(a, b)
                @test !lex(b, a)
                @test !dl(b, a)
            end

            2 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 1])
            b = monom_construct_from_vector(EV{T}, [0, 2])
            @test drl(b, a)
            @test !drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test !lex(a, b)
                @test !dl(a, b)
                @test lex(b, a)
                @test dl(b, a)
            end

            2 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 1])
            b = monom_construct_from_vector(EV{T}, [1, 1])
            @test !drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test !lex(a, b)
                @test !dl(a, b)
            end

            2 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 1])
            b = monom_construct_from_vector(EV{T}, [2, 2])
            @test !drl(b, a)
            @test drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test lex(a, b)
                @test dl(a, b)
                @test !lex(b, a)
                @test !dl(b, a)
            end

            3 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 0, 2])
            b = monom_construct_from_vector(EV{T}, [2, 0, 1])
            @test drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test lex(a, b)
                @test dl(a, b)
            end

            3 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, [1, 5, 0])
            b = monom_construct_from_vector(EV{T}, [1, 0, 1])
            @test !drl(a, b)

            25 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, ones(UInt, 25))
            b = monom_construct_from_vector(EV{T}, ones(UInt, 25))
            @test !drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test !lex(a, b)
                @test !dl(a, b)
            end

            30 > Groebner.monom_max_vars(EV{T}) && continue
            a = monom_construct_from_vector(EV{T}, ones(UInt, 30))
            b = monom_construct_from_vector(EV{T}, ones(UInt, 30))
            @test !drl(a, b)
            if !(EV{T} <: Groebner.AbstractPackedTuple{T})
                @test !lex(a, b)
                @test !dl(a, b)
            end
        end

        # test that different implementations agree
        for n in 1:10
            k = rand(1:60)

            implementations_to_test_local = filter(
                EV ->
                    !(EV{T} <: Groebner.AbstractPackedTuple{T}) &&
                        Groebner.monom_max_vars(EV{T}) >= k,
                implementations_to_test
            )

            t = div(typemax(UInt8), k) - 1
            x, y = rand(1:t, k), rand(1:t, k)
            if sum(x) >= Groebner.monom_overflow_threshold(UInt8)
                continue
            end
            if sum(y) >= Groebner.monom_overflow_threshold(UInt8)
                continue
            end
            as = [
                monom_construct_from_vector(EV{T}, x) for
                EV in implementations_to_test_local
            ]
            bs = [
                monom_construct_from_vector(EV{T}, y) for
                EV in implementations_to_test_local
            ]

            @test length(unique(map(lex, as, bs))) == 1
            @test length(unique(map(dl, as, bs))) == 1
            @test length(unique(map(drl, as, bs))) == 1

            # test that a < b && b < a does not happen
            for (a, b) in zip(as, bs)
                k > Groebner.monom_max_vars(a) && continue
                if lex(a, b)
                    @test !lex(b, a)
                end
                if dl(a, b)
                    @test !dl(b, a)
                end
                if drl(a, b)
                    @test !drl(b, a)
                end
            end
        end
    end
end

function test_circular_shift(a, b, n, Ord, answers)
    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:n]...]
    vars_to_index = Dict(x .=> 1:n)
    orders = map(
        i -> Groebner.ordering_transform(Ord(circshift(x, -i)), vars_to_index),
        0:(n - 1)
    )
    cmps = map(ord -> ((x, y) -> Groebner.monom_isless(x, y, ord)), orders)

    for (cmp, answer) in zip(cmps, answers)
        @test answer == cmp(a, b) && !answer == cmp(b, a)
    end
end

@testset "monoms, variable permutation" begin
    for T in (UInt64, UInt32, UInt16)
        for EV in implementations_to_test
            if EV{T} <: Groebner.AbstractPackedTuple
                continue
            end

            n = 3
            n >= Groebner.monom_max_vars(EV{T}) && continue

            a = monom_construct_from_vector(EV{T}, [5, 5, 3])
            b = monom_construct_from_vector(EV{T}, [1, 1, 10])

            test_circular_shift(a, b, n, Groebner.Lex, [false, false, true])
            test_circular_shift(a, b, n, Groebner.DegLex, [false, false, false])
            test_circular_shift(a, b, n, Groebner.DegRevLex, [false, false, false])

            n = 5
            n >= Groebner.monom_max_vars(EV{T}) && continue

            a = monom_construct_from_vector(EV{T}, [1, 2, 3, 4, 5])
            b = monom_construct_from_vector(EV{T}, [4, 3, 2, 1, 5])

            test_circular_shift(a, b, n, Groebner.Lex, [true, true, false, false, true])
            test_circular_shift(a, b, n, Groebner.DegLex, [true, true, false, false, true])
            test_circular_shift(
                a,
                b,
                n,
                Groebner.DegRevLex,
                [true, false, false, true, true]
            )
        end
    end
end

function test_orderings(n, v1, v2, ords_to_test)
    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:n]...]
    var_to_index = Dict(x .=> 1:n)
    for wo in ords_to_test
        ord = wo.ord
        ans = wo.ans
        internal_ord = Groebner.ordering_transform(ord, var_to_index)
        @test Groebner.monom_isless(v1, v2, internal_ord) == ans[1]
        @test Groebner.monom_isless(v2, v1, internal_ord) == ans[2]
        @test Groebner.monom_isless(v2, v2, internal_ord) == false
        @test Groebner.monom_isless(v1, v1, internal_ord) == false
    end
end

@testset "monom orders: WeightedOrdering" begin
    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:3]...]

    pv = Groebner.ExponentVector{T} where {T}

    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [3, 2, 1])

    @test_throws DomainError Groebner.WeightedOrdering(Dict(x .=> [-1, 0, 0]))

    ords_to_test = [
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 1, 1])), ans=[true, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 0, 1])), ans=[false, true]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 1, 0])), ans=[true, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 0, 0])), ans=[true, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 1, 5])), ans=[false, true]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 0, 0])), ans=[true, false])
    ]

    test_orderings(3, v1, v2, ords_to_test)

    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3])

    ords_to_test = [
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 1, 1])), ans=[false, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 0, 1])), ans=[false, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 1, 0])), ans=[false, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 0, 0])), ans=[false, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 1, 5])), ans=[false, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 0, 0])), ans=[false, false])
    ]

    test_orderings(3, v1, v2, ords_to_test)

    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:6]...]

    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3, 0, 0, 7])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [4, 5, 0, 0, 1, 4])

    ords_to_test = [
        (ord=Groebner.WeightedOrdering(Dict(x .=> [1, 1, 1, 1, 1, 1])), ans=[true, false]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 0, 0, 0, 0, 4])), ans=[false, true]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 2, 5, 0, 0, 0])), ans=[false, true]),
        (ord=Groebner.WeightedOrdering(Dict(x .=> [0, 2, 2, 0, 0, 0])), ans=[true, false])
    ]

    test_orderings(6, v1, v2, ords_to_test)
end

@testset "monom orders: ProductOrdering" begin
    pv = Groebner.ExponentVector{T} where {T}

    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:3]...]
    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [3, 2, 1])

    ords_to_test = [
        (
            ord=Groebner.ProductOrdering(Groebner.Lex(x[1]), Groebner.DegLex(x[2], x[3])),
            ans=[true, false]
        ),
        (
            ord=Groebner.ProductOrdering(Groebner.DegRevLex(x[1:2]), Groebner.Lex(x[3])),
            ans=[true, false]
        )
    ]

    test_orderings(3, v1, v2, ords_to_test)

    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:6]...]
    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [4, 1, 7, 0, 9, 8])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 6, 3, 5, 9, 100])

    ords_to_test = [
        (
            ord=Groebner.ProductOrdering(Groebner.DegLex(x[1]), Groebner.DegLex(x[2:6])),
            ans=[false, true]
        ),
        (
            ord=Groebner.ProductOrdering(Groebner.DegLex(x[1:2]), Groebner.DegLex(x[3:6])),
            ans=[true, false]
        ),
        (
            ord=Groebner.ProductOrdering(Groebner.DegLex(x[1:3]), Groebner.DegLex(x[4:6])),
            ans=[false, true]
        ),
        (
            ord=Groebner.ProductOrdering(Groebner.DegLex(x[1:4]), Groebner.DegLex(x[5:6])),
            ans=[true, false]
        ),
        (
            ord=Groebner.ProductOrdering(
                Groebner.ProductOrdering(Groebner.DegLex(x[1:2]), Groebner.DegLex(x[3:4])),
                Groebner.DegLex(x[5:6])
            ),
            ans=[true, false]
        ),
        (
            ord=Groebner.ProductOrdering(
                Groebner.ProductOrdering(Groebner.DegLex(x[1:3]), Groebner.DegLex(x[4])),
                Groebner.DegLex(x[5:6])
            ),
            ans=[false, true]
        ),
        (
            ord=Groebner.DegLex(x[1:3]) * Groebner.DegLex(x[4]) * Groebner.DegLex(x[5:6]),
            ans=[false, true]
        )
    ]

    test_orderings(6, v1, v2, ords_to_test)
end

@testset "monom orders: MatrixOrdering" begin
    R, x = AbstractAlgebra.QQ[["x$i" for i in 1:3]...]

    pv = Groebner.ExponentVector{T} where {T}

    v1 = Groebner.monom_construct_from_vector(pv{UInt64}, [1, 2, 3])
    v2 = Groebner.monom_construct_from_vector(pv{UInt64}, [3, 2, 1])

    ord1 = Groebner.MatrixOrdering(
        x,
        [
            1 0 0
            0 1 0
            0 0 1
        ]
    )
    ord2 = Groebner.MatrixOrdering(
        x,
        [
            1 0 2;
        ]
    )
    ord3 = Groebner.MatrixOrdering(
        x,
        [
            0 0 0
            0 1 0
            1 1 1
        ]
    )
    @test_throws DomainError Groebner.MatrixOrdering(
        x,
        [
            1 0;
        ]
    )

    mo_to_test = [
        (ord=ord1, ans=[true, false]),
        (ord=ord2, ans=[false, true]),
        (ord=ord3, ans=[false, false])
    ]
    test_orderings(3, v1, v2, mo_to_test)
end
