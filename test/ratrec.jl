using Test, Groebner

@testset "CRT table" begin
    table_qq = [Rational{BigInt}[], Rational{BigInt}[0,0,0], Rational{BigInt}[0,0]]
    res = [Rational{BigInt}[], Rational{BigInt}[2//3, 0//1, -5//8], Rational{BigInt}[-11//1, 7//4]]
    
    modulo = BigInt(131)
    table_zz = [BigInt[mod(numerator(r) * invmod(denominator(r), modulo), modulo) for r in row] for row in res]
    mask = [falses(length(table_qq[i])) for i in 1:length(table_qq)]
    success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo, mask)
    @test !success

    modulo = BigInt(251)
    table_zz = [BigInt[mod(numerator(r) * invmod(denominator(r), modulo), modulo) for r in row] for row in res]
    success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo, mask)
    @test success && table_qq == res
    for tasks in [1,2,3,4]
        success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo, mask, tasks=tasks)
        @test success && table_qq == res
    end
end
