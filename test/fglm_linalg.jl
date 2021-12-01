

using .GroebnerBases: findpivot, rowreduce!, linear_relation!

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

@testset "fglm linear recurrence" begin
    @test findpivot([[1]], 1, 1) == (1, 1)
    @test findpivot([[0]], 1, 1) == (0, 0)

    A = [
        QQ.([0, 1]),
        QQ.([1, 0])
    ]
    @test findpivot(A, 1, 1) == (2, 1)
    @test findpivot(A, 1, 2) == (1, 2)

    A = [
        [1, 0, 0],
        [0, 0, 2]
    ]
    @test findpivot(A, 2, 2) == (2, 3)

    #=
    v = [2, 2]
    Ar, vector, ps = GroebnerBases.rowreduce!(A, v)
    @test ps == [(1, 1), (2, 3)] && A == [[1, 0, 0], [0, 0, 1]]

    A = [
        QQ.([0, 0, 8, 0]),
        QQ.([0, 0, 0, 5]),
        QQ.([0, 1, 0, 0])
    ]

    v = QQ.([0, 1, 1, 1])
    Ar, v, ps = GroebnerBases.rowreduce!(A, v)
    @test ps == [(1, 1), (2, 2), (3, 3)] && Ar == [
        QQ.([0, 1, 0, 0]),
        QQ.([0, 0, 1, 0]),
        QQ.([0, 0, 0, 1])
    ] && v == [1//8, 1//5, ]

    @test roworder == [3, 4, 2, 1]

    A = [
        QQ.([0, 0, 12, 0]),
        QQ.([-1, 1, 11, 0]),
        QQ.([0, 0, 6, 0]),
    ]
    Ar, ps, roworder = GroebnerBases.rowreduce!(A)
    @test ps == [(1, 1), (2, 3)] && Ar == [
        QQ.([-1,1, 0, 0]),
        QQ.([0, 0, 12, 0]),
        QQ.([0, 0, 0, 0]),
    ]
    @test roworder == [2, 1, 3]

    A = [
        QQ.([0, 0, 1]),
        QQ.([0, 2, 0]),
        QQ.([-3, 0, 0])
    ]
    Ar, ps, roworder = GroebnerBases.rowreduce!(A)
    @test ps == [(1, 1), (2, 2), (3, 3)] && roworder == [3, 2, 1]

    =#

    A = [
        QQ.([0, 0, 0]),
        QQ.([0, 0, 1]),
        QQ.([8, 0, 0]),
        QQ.([0, 5, 0])
    ]
    v = QQ.([1, -2, 3])
    λ, ex = GroebnerBases.linear_relation!(deepcopy(A), v)
    println(λ)
    @test λ == QQ.([0, 3, 1//8, -2//5]) && ex


    A = [
        QQ.([0, 1, 0]),
        QQ.([8, 0, 4]),
        QQ.([0, 0, 0]),
        QQ.([0, 5, 0])
    ]
    v = QQ.([-2, 0, -1])
    λ, ex = GroebnerBases.linear_relation!(deepcopy(A), v)
    @test λ == QQ.([0, -1//4, 0, 0]) && ex


    A = [
        QQ.([1, 3]),
        QQ.([2, 4]),
    ]
    v = QQ.([4, 10])
    λ, ex = GroebnerBases.linear_relation!(A, v)
    @test λ == QQ.([2, 1]) && ex

    A = [
        QQ.([0, 0, 0, 0, 2]),
        QQ.([0, 0, 1, 0, 0]),
        QQ.([0, 0, 1, 0, 0]),
        QQ.([1, 2, 0, 1, 0]),
    ]
    v = QQ.([-2, -4, 1, -2, 8])
    λ, ex = GroebnerBases.linear_relation!(A, v)
    @test λ == QQ.([4, 1, 0, -2])

    A = [[0]]
    v = [1]
    λ, ex = GroebnerBases.linear_relation!(A, v)
    @test !ex

    A = [QQ.([1, 2]), QQ.([2, 4])]
    v = QQ.([0, 1])
    λ, ex = GroebnerBases.linear_relation!(A, v)
    @test !ex
end
