using AbstractAlgebra
using Combinatorics
import Random

#=
TODO: curious case:

system = AbstractAlgebra.Generic.MPoly{AbstractAlgebra.GFElem{Int64}}[
    x^2 + 1, 
    x*y, 
    y*z + 1
] 
ord = Groebner.DegLex{Nothing}(nothing, false))

y*z + 1
x*y
x^2 + 1
=#

@testset "homogenization, simple" failfast = true begin
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

        # Test that the basis obtained with the use of homogenization
        # *coincides* with the one obtained without it
        R, (x, y, z) = PolynomialRing(field, ["x", "y", "z"])
        for case in [
            (system=[x, y, z], ord=Groebner.Lex()),
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
            (system=Groebner.katsuran(4, ground=field), ord=Groebner.DegRevLex()),
            (system=Groebner.rootn(5, ground=field), ord=Groebner.Lex()),
            (system=Groebner.sparse5(ground=field), ord=Groebner.Lex()),
            (system=Groebner.reimern(5, ground=field), ord=Groebner.DegRevLex())
        ]
            gb1 = Groebner.groebner(case.system, ordering=case.ord, homogenize=:no)
            gb2 = Groebner.groebner(case.system, ordering=case.ord, homogenize=:yes)
            @info "" gb1 gb2 case field
            # @test Groebner.isgroebner(gb1, ordering=case.ord)
            @test gb1 == gb2
        end
    end
end
