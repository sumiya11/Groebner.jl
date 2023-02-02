using AbstractAlgebra
using Base.Threads

@testset "multi-threading" begin
    @warn "Using $(nthreads()) threads."
    for ground in [QQ, GF(2^31-1)]
        R, (x,y,z) = PolynomialRing(ground, ["x","y","z"], ordering=:degrevlex)
        cases = [
            [R(1)],
            [x, y, z, y, y, z, x, x, x, x, y],
            [(x + y + z)^5,
            (x + y)^4,
            (x + z)^7],
            [x^2*y + BigInt(10)^100, x*y^2 + BigInt(11)^100],
            [x^2*y + BigInt(10)^500, x*y^2 + BigInt(11)^500],
            [x^2*y + BigInt(10)^1000, x*y^2 + BigInt(11)^1000],
            Groebner.cyclicn(5, ground=ground),
            Groebner.cyclicn(6, ground=ground),
            Groebner.katsuran(6, ground=ground),
            Groebner.katsuran(7, ground=ground)
        ]
        for case in cases
            gb1 = Groebner.groebner(case, threading=false, ordering=:degrevlex)
            gb2 = Groebner.groebner(case, threading=true, ordering=:degrevlex)
            @test gb1 == gb2
        end
    end
end
