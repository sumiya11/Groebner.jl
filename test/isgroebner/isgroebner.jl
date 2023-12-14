using AbstractAlgebra

@testset "isgroebner" begin
    R, x = polynomial_ring(GF(2^31 - 1), "x")
    @test Groebner.isgroebner([x])
    @test Groebner.isgroebner([x, x, x, x])
    @test !Groebner.isgroebner([x^2, x^2 + 1])

    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)

    @test Groebner.isgroebner([R(1)])
    @test Groebner.isgroebner([x])
    @test Groebner.isgroebner([x, y])
    @test Groebner.isgroebner([x, x])
    @test Groebner.isgroebner([x, y, z])

    @test !Groebner.isgroebner([x + y, x])
    @test Groebner.isgroebner([z * x, z * x, R(1)])

    fs = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))

    fs = Groebner.change_ordering(Groebner.rootn(5, ground=GF(2^31 - 1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test Groebner.isgroebner(Groebner.groebner(fs, reduced=false))

    fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31 - 1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))

    fs = Groebner.change_ordering(Groebner.noonn(3, ground=GF(2^31 - 1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test Groebner.isgroebner(Groebner.groebner(fs, reduced=false))

    fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31 - 1)), :degrevlex)
    @test !Groebner.isgroebner(fs)
    @test Groebner.isgroebner(Groebner.groebner(fs))
end

@testset "isgroebner orderings" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])

    @test Groebner.isgroebner([x^2 + y, y])
    @test Groebner.isgroebner([x^2 + y, y], ordering=Groebner.Lex(x, y, z))
    @test !Groebner.isgroebner([x^2 + y, y], ordering=Groebner.Lex(y, x, z))

    @test Groebner.isgroebner([x^2 + y, y], ordering=Groebner.DegLex())
    @test Groebner.isgroebner([x^2 + y, y], ordering=Groebner.DegLex(y, x, z))
end

@testset "isgroebner certify" begin
    for certify in [false, true]
        for field in [GF(17), GF(2^31 - 1), QQ]
            R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], ordering=:degrevlex)

            @test Groebner.isgroebner([R(0)], certify=certify)
            @test Groebner.isgroebner([R(1)], certify=certify)
            @test Groebner.isgroebner([x], certify=certify)
            @test Groebner.isgroebner([x, y], certify=certify)
            @test Groebner.isgroebner([x, x], certify=certify)
            @test Groebner.isgroebner([x, y, z], certify=certify)

            @test !Groebner.isgroebner([x + y, x], certify=certify)
            @test Groebner.isgroebner([z * x, z * x, R(1)], certify=certify)
        end
        fs = Groebner.rootn(3, ground=GF(2^31 - 1), ordering=:degrevlex)
        @test !Groebner.isgroebner(fs, certify=certify)
        @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

        fs = Groebner.rootn(3, ground=QQ, ordering=:degrevlex)
        @test !Groebner.isgroebner(fs, certify=certify)
        @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

        fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31 - 1)), :degrevlex)
        @test !Groebner.isgroebner(fs, certify=certify)
        @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

        # TODO: figure out why this is broken!
        # fs = Groebner.change_ordering(Groebner.noonn(2, ground=QQ), :degrevlex)
        # @test !Groebner.isgroebner(fs, certify=certify)
        # @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

        fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31 - 1)), :degrevlex)
        @test !Groebner.isgroebner(fs, certify=certify)
        @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)
    end
end
