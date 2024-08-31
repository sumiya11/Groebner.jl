using Base.Threads

@testset "groebner, multi-threading, Zp" begin
    @info "Testing multi-threading over Zp using $(nthreads()) threads"

    R, (x, y, z) = GF(2^31 - 1)["x", "y", "z"]

    s = [R(1), R(2), R(4)]
    @test Groebner.groebner(s, threaded=:yes) ==
            Groebner.groebner(s, threaded=:no) ==
            Groebner.groebner(s, threaded=:auto) ==
            [R(1)]

    s = [x^(2^10) + 1, y^3000 + 2, z^1000 + 3]
    @test Groebner.groebner(s, threaded=:yes) ==
            Groebner.groebner(s, threaded=:no) ==
            Groebner.groebner(s, threaded=:auto) ==
            [z^1000 + 3, y^3000 + 2, x^(2^10) + 1]

    linalg = [:deterministic, :randomized]
    threaded = [:yes, :no, :auto]
    grid = [(linalg=l, threaded=t) for l in linalg for t in threaded]
    for system in [
        Groebner.Examples.katsuran(3, internal_ordering=:lex, k=GF(2^31 - 1)),
        Groebner.Examples.katsuran(4, internal_ordering=:lex, k=GF(2^20 + 7)),
        Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
        Groebner.Examples.eco5(internal_ordering=:deglex, k=GF(2^20 + 7)),
        Groebner.Examples.eco5(internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.ojika4(internal_ordering=:lex, k=GF(2^20 + 7)),
        Groebner.Examples.henrion5(internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=GF(2^10 + 7)),
        Groebner.Examples.ku10(internal_ordering=:degrevlex, k=GF(2^10 + 7))
    ]
        results = Array{Any}(undef, size(grid))
        for (i, kw) in enumerate(grid)
            gb = Groebner.groebner(system; kw...)
            results[i] = gb
        end
        @test length(unique(results)) == 1
        @test Groebner.isgroebner(results[1])
    end
end

@testset "groebner, multi-threading, QQ" begin
    @info "Testing multi-threading over QQ using $(nthreads()) threads"

    R, (x, y) = QQ["x", "y"]

    threaded = [:yes, :no, :auto]
    grid = [(threaded=t,) for t in threaded]
    for system in [
        [x - 1, y - 2],
        [x + (BigInt(2)^1000 + 1) // 2^61, x * y + BigInt(2)^(2^10)],
        Groebner.Examples.katsuran(3, internal_ordering=:lex, k=QQ),
        Groebner.Examples.katsuran(4, internal_ordering=:lex, k=QQ),
        Groebner.Examples.eco5(internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.ojika4(internal_ordering=:lex, k=QQ),
        Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.ku10(internal_ordering=:degrevlex, k=QQ),
        [x + BigInt(2)^1_000 // (BigInt(2)^1_000 - 1), y - BigInt(2)^1_000],
        [x + BigInt(2)^10_000 // (BigInt(2)^10_000 - 1), y^2 - BigInt(2)^10_000 * y]
    ]
        results = Array{Any}(undef, size(grid))
        for (i, kw) in enumerate(grid)
            gb = Groebner.groebner(system; kw...)
            results[i] = gb
        end
        @test length(unique(results)) == 1
        @test Groebner.isgroebner(results[1])
    end
end
