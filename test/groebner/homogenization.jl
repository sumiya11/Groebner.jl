using AbstractAlgebra
using Combinatorics
import Random

@testset "homogenization, basic" begin
    for field in [GF(113), GF(2^62 + 135), QQ]
        R, x = polynomial_ring(field, "x")
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:yes) == [x + 1]
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:no) == [x + 1]
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:auto) == [x + 1]

        for ordering in [:lex, :deglex, :degrevlex]
            R, (x, y) = polynomial_ring(field, ["x", "y"], ordering=ordering)

            @test Groebner.groebner([R(3)], homogenize=:yes) == [R(1)]
            @test Groebner.groebner([R(0)], homogenize=:yes) == [R(0)]
            @test Groebner.groebner([x], homogenize=:yes) == [x]
            @test Groebner.groebner([x + 11y], homogenize=:yes) == [x + 11y]
            @test Groebner.groebner([x, y], homogenize=:yes) == [y, x]
            @test Groebner.isgroebner(Groebner.groebner([x + 1, y + 2], homogenize=:yes))

            if ordering === :lex
                gb = Groebner.groebner([x + y^2, x^2 - 1], homogenize=:yes)
                @test (y^4 - 1) in gb && (x + y^2) in gb
            end

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

@testset "homogenization, orderings" begin
    for field in [GF(113), GF(2^62 + 135), QQ]
        # Test that the basis obtained with the use of homogenization
        # *coincides* with the one obtained without it
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"])
        for case in [
            (system=[R(1)], ord=Groebner.Lex()),
            (system=[R(5)], ord=Groebner.DegRevLex()),
            (system=[x, y, z], ord=Groebner.Lex()),
            (system=[x, y, z], ord=Groebner.Lex(z) * Groebner.Lex(x, y)),
            (system=[x^2 + 1, x * y + 2, y * z + 3], ord=Groebner.DegLex()),
            (
                system=[x^20 - x^15 - x^5 + y * z, x * y^10 + x^3 * y^3 + x * y],
                ord=Groebner.DegRevLex()
            ),
            (system=[x + 5y, x + 7z, y + 11z], ord=Groebner.Lex(z) * Groebner.Lex(x, y)),
            (system=[x + 5y + 11z], ord=Groebner.DegRevLex(y) * Groebner.Lex(z, x)),
            (system=[x + 5y + 11z], ord=Groebner.DegRevLex(z) * Groebner.DegRevLex(y, x)),
            (system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1], ord=Groebner.Lex()),
            (
                system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1],
                ord=Groebner.DegRevLex()
            ),
            (
                system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1],
                ord=Groebner.DegRevLex(z) * Groebner.DegRevLex(x, y)
            ),
            (
                system=[x * y^2 + x + 1, y * z^2 + y + 1, x * y * z^4 - x * z^2 - y + z],
                ord=Groebner.DegRevLex(y) * Groebner.DegRevLex(z, x)
            ),
            (system=Groebner.katsuran(4, ground=field), ord=Groebner.DegRevLex()),
            (system=Groebner.rootn(5, ground=field), ord=Groebner.Lex()),
            (system=Groebner.sparse5(ground=field), ord=Groebner.Lex()),
            (system=Groebner.reimern(5, ground=field), ord=Groebner.DegRevLex())
        ]
            ord = case.ord
            system = case.system
            @debug "" system ord field

            gb1 = Groebner.groebner(system, ordering=ord, homogenize=:no)
            gb2 = Groebner.groebner(system, ordering=ord, homogenize=:yes)

            @test Groebner.isgroebner(gb1, ordering=ord)

            @test gb1 == gb2

            # Also test learn / apply
            if field != QQ
                context3, gb3 =
                    Groebner.groebner_learn(system, ordering=ord, homogenize=:yes)
                context4, gb4 =
                    Groebner.groebner_learn(system, ordering=ord, homogenize=:no)
                for _ in 1:4
                    flag3, gb33 = Groebner.groebner_apply!(context3, system)
                    flag4, gb44 = Groebner.groebner_apply!(context4, system)
                    @test flag3 && flag4
                    @test gb1 == gb3 == gb4 == gb33 == gb44
                end
            end
        end
    end
end
