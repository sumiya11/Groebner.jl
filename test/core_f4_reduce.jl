
@testset "ff f4 reduce sort invariant" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x","y","z"], ordering=:lex)

    @test Groebner.groebner([x, y]) == Groebner.groebner([y, x]) == [y, x]

    fs = [x^2 + y, x*y]
    @test Groebner.groebner(fs) == [y^2, x*y, x^2 + y]

end

@testset "ff f4 yesreduce" begin

    root = Groebner.change_ordering(Groebner.rootn(3, ground=GF(2^31 - 1)), :degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, reduced=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2*x3 + x3^2, x3^3 + 2147483646]

    root = Groebner.change_ordering(Groebner.rootn(6, ground=GF(2^31 - 1)), :degrevlex)
    x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
    gb = Groebner.groebner(root, reduced=true)
    @test gb == [
        x1 + x2 + x3 + x4 + x5 + x6,
        x2^2 + x2*x3 + x3^2 + x2*x4 + x3*x4 + x4^2 +
        x2*x5 + x3*x5 + x4*x5 + x5^2 + x2*x6 + x3*x6 + x4*x6 +
        x5*x6 + x6^2, x3^3 + x3^2*x4 + x3*x4^2 + x4^3 +
        x3^2*x5 + x3*x4*x5 + x4^2*x5 + x3*x5^2 + x4*x5^2 +
        x5^3 + x3^2*x6 + x3*x4*x6 + x4^2*x6 + x3*x5*x6 +
        x4*x5*x6 + x5^2*x6 + x3*x6^2 + x4*x6^2 + x5*x6^2 +
        x6^3, x4^4 + x4^3*x5 + x4^2*x5^2 + x4*x5^3 + x5^4 +
        x4^3*x6 + x4^2*x5*x6 + x4*x5^2*x6 + x5^3*x6 + x4^2*x6^2 +
        x4*x5*x6^2 + x5^2*x6^2 + x4*x6^3 + x5*x6^3 + x6^4, x5^5 +
        x5^4*x6 + x5^3*x6^2 + x5^2*x6^3 + x5*x6^4 + x6^5,
        x6^6 + 2147483646]

    ku = Groebner.change_ordering(Groebner.ku10(ground=GF(2^31-1)), :degrevlex)
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = gens(parent(first(ku)))
    gb = Groebner.groebner(ku, reduced=true)

    @test gb == [
        x9 + 1272065637*x10 + 875418006,
        x8 + 1529540685*x10 + 617942964,
        x7 + 1539832471*x10 + 607651173,
        x6 + 1314302432*x10 + 833181218,
        x5 + 1453635454*x10 + 693848197,
        x4 + 673118236*x10 + 1474365406,
        x3 + 269783061*x10 + 1877700587,
        x2 + 1042807874*x10 + 1104675778,
        x1 + 389079675*x10 + 1758403970,
        x10^2 + 1222705397*x10 + 924778249]

end
