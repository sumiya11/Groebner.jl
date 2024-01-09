using AbstractAlgebra

@testset "fglm, Zp" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"])

    @test Groebner.fglm([one(R)], Groebner.DegRevLex(), Groebner.Lex()) ==
          Groebner.fglm([one(R)], Groebner.DegLex(), Groebner.Lex()) ==
          [one(R)]

    sys = [x, y]
    @test Groebner.fglm(sys, Groebner.DegRevLex(), Groebner.Lex()) ==
          Groebner.fglm(sys, Groebner.DegLex(), Groebner.Lex()) ==
          [y, x]

    noon = [10x * y^2 - 11 * x + 10, 10x^2 * y - 11y + 10]

    gb = Groebner.groebner(noon, ordering=Groebner.DegRevLex())
    true_gb_lex = [
        y^5 + 1952257860 * y^4 + 1288490186 * y^3 + 2 * y^2 + 1028839893 * y + 644245093,
        x + y^4 + 1952257860 * y^3 + 644245093 * y^2 + y + 1952257860
    ]
    gb_fglm = Groebner.fglm(gb, Groebner.DegRevLex(), Groebner.Lex())
    @test gb_fglm == true_gb_lex
    @test parent(first(gb_fglm)) == R

    gb = Groebner.groebner(noon, ordering=Groebner.DegLex())
    @test Groebner.fglm(gb, Groebner.DegRevLex(), Groebner.Lex()) == true_gb_lex

    gb = Groebner.groebner(noon, ordering=Groebner.Lex())
    @test Groebner.fglm(gb, Groebner.Lex(), Groebner.Lex()) == true_gb_lex
end

@testset "fglm, the rationals" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering=:degrevlex)

    gb = [
        z^2 - (9 // 490) * y - 201 // 980 * z + 13 // 980,
        y * z - 4 // 35 * y + 2 // 35 * z - 1 // 35,
        y^2 - 1 // 10 * y - 6 // 5 * z + 1 // 10,
        x^2 - 4 // 5 * y + 2 // 5 * z - 1 // 5
    ]

    gb_fglm = Groebner.fglm(gb, Groebner.DegRevLex(), Groebner.Lex())
    @test parent(first(gb_fglm)) == R
    @test Groebner.fglm(gb, Groebner.DegRevLex(), Groebner.Lex()) == [
        z^3 - 313 // 980 * z^2 + 37 // 980 * z - 1 // 490,
        y - 490 // 9 * z^2 + 67 // 6 * z - 13 // 18,
        x^2 - 392 // 9 * z^2 + 28 // 3 * z - 7 // 9
    ]
end

@testset "fglm, non-shape" begin
    R, (x, y, z, t) =
        AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x", "y", "z", "t"])
    sys = [
        y^2 * z + 2 * x * y * t - 2 * x - z,
        -x^3 * z + 4 * x * y^2 * z + 4 * x^2 * y * t + 2 * y^3 * t + 4 * x^2 - 10 * y^2 + 4 * x * z - 10 * y * t + 2,
        2 * y * z * t + x * t^2 - x - 2 * z,
        -x * z^3 + 4 * y * z^2 * t + 4 * x * z * t^2 + 2 * y * t^3 + 4 * x * z + 4 * z^2 - 10 * y * t - 10 * t^2 + 2
    ]

    gb_lex = Groebner.groebner(sys, ordering=Groebner.Lex())
    gb_drl = Groebner.groebner(sys, ordering=Groebner.DegRevLex())
    gb_fglm = Groebner.fglm(gb_drl, Groebner.DegRevLex(), Groebner.Lex())

    #! format: off
    @test gb_fglm == gb_lex == [
        t^20 - 191//3*t^16 - 8000//27*t^14 - 11746//27*t^12 + 11746//27*t^8 + 8000//27*t^6 + 191//3*t^4 - 1,
        z*t^14 + 9*z*t^12 + 79//3*z*t^10 + 559//27*z*t^8 - 559//27*z*t^6 - 79//3*z*t^4 - 9*z*t^2 - z,
        z^2*t^6 + 7//3*z^2*t^4 - 7//3*z^2*t^2 - z^2 - 9//320*t^18 + 171//3200*t^16 + 573//320*t^14 + 23413//4800*t^12 - 3427//960*t^10 - 1861//96*t^8 - 1219//960*t^6 + 71587//4800*t^4 + 1477//480*t^2 - 1471//3200,
        z^4*t^2 - z^4 + 4*z^2*t^4 - 4*z^2 + 2943//16000*t^18 + 297//16000*t^16 - 11709//1000*t^14 - 111367//2000*t^12 - 137033//1600*t^10 - 13843//1600*t^8 + 166293//2000*t^6 + 8312//125*t^4 + 224387//16000*t^2 - 34867//16000,
        z^9 - 16*z^7 + 208//3*z^3*t^4 + 224//3*z^3*t^2 + 112*z^3 - 261//8*z*t^12 - 1125//4*z*t^10 - 6179//8*z*t^8 - 2785//6*z*t^6 + 14735//24*z*t^4 + 6049//12*z*t^2 + 1411//8*z,
        y - 1//128*z^8*t + 1//8*z^6*t - 13//24*z^2*t^5 - 7//12*z^2*t^3 - 7//8*z^2*t - 41841//409600*t^19 - 63477//2048000*t^17 + 668091//102400*t^15 + 49508897//1536000*t^13 + 32128219//614400*t^11 + 973303//122880*t^9 - 15595723//307200*t^7 - 60484147//1536000*t^5 - 1953023//245760*t^3 + 523427//2048000*t,
        x - 1//96*z^7 + 7//48*z^5 - 5//32*z^3*t^4 - 5//48*z^3*t^2 + 41//96*z^3 - 693//4096*z*t^12 - 2889//2048*z*t^10 - 14955//4096*z*t^8 - 2415//1024*z*t^6 + 8197//4096*z*t^4 + 14821//6144*z*t^2 + 7547//4096*z
    ]
    #! format: on
end
