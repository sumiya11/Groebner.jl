using Base.Threads

@testset "groebner multi-threading Zp" begin
    if Groebner._threaded[]
        @info "Testing multi-threading with $(nthreads()) threads"

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

        linalgs = [:deterministic, :randomized]
        threadeds = [:yes, :no, :auto]
        grid = [(linalg=l, threaded=t) for l in linalgs for t in threadeds]
        for system in [
            Groebner.katsuran(3, ordering=:lex, ground=GF(2^31 - 1)),
            Groebner.katsuran(4, ordering=:lex, ground=GF(2^20 + 7)),
            Groebner.eco5(ordering=:deglex, ground=GF(2^20 + 7)),
            Groebner.eco5(ordering=:degrevlex, ground=GF(2^20 + 7)),
            Groebner.cyclicn(5, ordering=:degrevlex, ground=GF(2^20 + 7)),
            Groebner.cyclicn(6, ordering=:degrevlex, ground=GF(2^20 + 7)),
            Groebner.ojika4(ordering=:lex, ground=GF(2^20 + 7)),
            Groebner.henrion5(ordering=:degrevlex, ground=GF(2^20 + 7)),
            Groebner.noonn(6, ordering=:degrevlex, ground=GF(2^10 + 7)),
            Groebner.ku10(ordering=:degrevlex, ground=GF(2^10 + 7)),
            Groebner.sparse5(ordering=:degrevlex, ground=GF(2^10 + 7))
        ]
            results = Array{Any}(undef, size(grid))
            for (i, kw) in enumerate(grid)
                gb = Groebner.groebner(system; kw...)
                results[i] = gb
            end
            @test length(unique(results)) == 1
            @test Groebner.isgroebner(results[1])
        end
    else
        @info "Skipping multi-threading tests"

        R, (x, y, z) = GF(2^31 - 1)["x", "y", "z"]
        @test_throws AssertionError Groebner.groebner([x], threaded=:yes)
    end
end
