using AbstractAlgebra
using Combinatorics
import Random

@testset "homogenization simple" failfast = true begin
    for field in [GF(2), GF(2^62 + 135), QQ]
        for ordering in [:lex, :deglex, :degrevlex]
            R, (x, y) = PolynomialRing(field, ["x", "y"], ordering=ordering)

            @test Groebner.groebner([R(3)], homogenize=:yes) == [R(1)]
            @test Groebner.groebner([R(0)], homogenize=:yes) == [R(0)]
            @test Groebner.groebner([x], homogenize=:yes) == [x]
            @test Groebner.groebner([x + 11y], homogenize=:yes) == [x + 11y]
            @test Groebner.groebner([x, y], homogenize=:yes) == [y, x]
            @test Groebner.isgroebner(Groebner.groebner([x + 1, y + 2], homogenize=:yes))

            for case in [
                Groebner.katsuran(2, ground=field, ordering=ordering),
                Groebner.katsuran(3, ground=field, ordering=ordering),
                Groebner.katsuran(4, ground=field, ordering=ordering),
                Groebner.noonn(2, ground=field, ordering=ordering),
                Groebner.noonn(3, ground=field, ordering=ordering)
            ]
                gb = Groebner.groebner(case, homogenize=:yes)
                @test Groebner.isgroebner(gb)
            end
        end
    end
end
