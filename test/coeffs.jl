
import Primes

@testset "Constants" begin
    @test Groebner.FIRST_COMPUTE_PRIME == 2^31-1
    @test Groebner.FIRST_CHECK_PRIME == 2^30 + 3
end

@testset "Lucky primes" begin
    v = [[BigInt(1), BigInt(2), BigInt(7)*123, BigInt(9)*(2^31-1)], [BigInt(11)*4]]
    vc = deepcopy(v)

    @test Groebner.isluckyprime(v, 5)
    @test !Groebner.isluckyprime(v, 11)
    @test !Groebner.isluckyprime(v, 2)
    @test !Groebner.isluckyprime(v, 7)
    @test !Groebner.isluckyprime(v, 2^31-1)
    @test Groebner.isluckyprime(v, 13)

    @test Groebner.nextluckyprime([BigInt[]]) == Groebner.FIRST_COMPUTE_PRIME
    @test Groebner.nextluckyprime(v) == Groebner.nextluckyprime(v, 2^31-1) == Primes.nextprime(2^31)
    @test Groebner.nextluckyprime([[Primes.nextprime(BigInt(2^31))]]) == 2^31-1
    @test Groebner.nextluckyprime([[Primes.nextprime(BigInt(2^31))]], 2^31-1) == Primes.nextprime(1+Primes.nextprime(BigInt(2^31)))

    @test v == vc
end

@testset "Good primes" begin
    v = [[BigInt(1), BigInt(2), BigInt(7)*123, BigInt(9)*(2^31-1)], [BigInt(11)*4]]
    vc = deepcopy(v)

    moduli = [17]
    @test Groebner.nextgoodprime(v, moduli, 16) == 19
    @test Groebner.nextgoodprime(v, [17, 19], 16) == 23
    @test Groebner.nextgoodprime(v, Int[], 17) == 19
    @test Groebner.nextgoodprime(v, Int[]) == Groebner.FIRST_CHECK_PRIME

    @test v == vc
end

@testset "Scaling coefficients" begin
    v = [
            [BigInt(1)//5, BigInt(2)//3, BigInt(7)//11],
            [BigInt(11), 1],
            [BigInt(-2)//1, BigInt(-1)//2],
            [BigInt(-1)//7, BigInt(3)//14, BigInt(-12)//7],
            [BigInt(1)//2],
            [BigInt(4115)//4115226, 1]
    ]
    vc = deepcopy(v)
    vv = Groebner.scale_denominators(v)

    @test vv == [
            [BigInt(33), BigInt(55)*2, BigInt(15)*7],
            [BigInt(11), 1],
            [BigInt(-4), BigInt(-1)],
            [BigInt(-2), BigInt(3), BigInt(-24)],
            [BigInt(1)],
            [BigInt(4115), 4115226]
    ]
    @test vc == v
end

@testset "Reconstructing coefficients modulo" begin
    v = [[BigInt(3), BigInt(4)], [BigInt(5), BigInt(6)]]
    qq = Vector{Vector{Rational{BigInt}}}(undef, 0)
    modulo = BigInt(11)
    vc = deepcopy(v)

    Groebner.reconstruct_modulo!(qq, v, modulo)
    @test v == vc
    @test qq == [[-2//3, 1//3], [-1//2, 1//2]]

    Groebner.reconstruct_modulo!(qq, v, modulo)
    @test v == vc
    @test qq == [[-2//3, 1//3], [-1//2, 1//2]]

    v = [[BigInt(9), BigInt(10)], [BigInt(11), BigInt(12)]]
    modulo = BigInt(37)
    vc = deepcopy(v)
    Groebner.reconstruct_modulo!(qq, v, modulo)
    @test v == vc
    @test qq == [[-1//4, 3//4], [-4//3, -1//3]]

end

@testset "Reconstructing coefficients crt" begin
    modulo = BigInt(1)
    ring_ff = Groebner.PolyRing(2, 3, :deglex, 5, :input)
    coeffs_ff = [[UInt64(1), UInt64(2)], [UInt64(4), UInt64(3)]]
    gbcoeffs_accum = Vector{Vector{BigInt}}(undef, 0)
    coeffs_ff_c = deepcopy(coeffs_ff)

    Groebner.reconstruct_crt!(gbcoeffs_accum, modulo, coeffs_ff, ring_ff)

    @test coeffs_ff_c == coeffs_ff
    @test gbcoeffs_accum == [[BigInt(1), 2], [BigInt(4), 3]]
    @test modulo == 5

    ring_ff = Groebner.PolyRing(2, 3, :deglex, 7, :input)
    coeffs_ff = [[UInt64(6), UInt64(1)], [UInt64(5), UInt64(3)]]
    coeffs_ff_c = deepcopy(coeffs_ff)

    Groebner.reconstruct_crt!(gbcoeffs_accum, modulo, coeffs_ff, ring_ff)
    @test coeffs_ff_c == coeffs_ff
    @test gbcoeffs_accum == [[BigInt(6), 22], [BigInt(19), 3]]
    @test modulo == 35
end
