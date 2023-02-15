
using AbstractAlgebra

@testset "isgroebner certify" begin
    for field in [GF(17), GF(2^31-1), QQ]
        for certify in [false, true]
            R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:degrevlex)

            @test Groebner.isgroebner([R(1)], certify=certify)
            @test Groebner.isgroebner([x], certify=certify)
            @test Groebner.isgroebner([x, y], certify=certify)
            @test Groebner.isgroebner([x, x], certify=certify)
            @test Groebner.isgroebner([x, y, z], certify=certify)

            @test !Groebner.isgroebner([x + y, x], certify=certify)
            @test Groebner.isgroebner([z*x, z*x, R(1)], certify=certify)

            fs = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31-1)), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

            fs = Groebner.change_ordering(Groebner.rootn(3, ground=QQ), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

            fs = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31-1)), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

            fs = Groebner.change_ordering(Groebner.noonn(2, ground=QQ), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

            fs = Groebner.change_ordering(Groebner.noonn(6, ground=GF(2^31-1)), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)

            fs = Groebner.change_ordering(Groebner.noonn(6, ground=QQ), :degrevlex)
            @test !Groebner.isgroebner(fs, certify=certify)
            @test Groebner.isgroebner(Groebner.groebner(fs), certify=certify)
        
        end
    end
end
