
using .GroebnerBases: reduce_modulo, scale_denominators,
                    rational_reconstruction, crt

using .GroebnerBases: Primes

R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"])
Rz, (xz, yz, zz, wz) = PolynomialRing(ZZ, ["x", "y", "z", "w"])


@testset "Scaling to integers" begin
    fs = [
        (1//2)x, (1//3)y, (1//4)z, (1//5)w
    ]
    scaled = scale_denominators(fs)
    Rzz = parent(first(scaled))
    @test scaled == gens(Rzz) && base_ring(Rzz) == ZZ

    fs = [
        (1//2)x^2 + (3//2)y^2 + (5//4)z,
        (3//2)x^2 + (5//3)y^2 + (7//5)z^2 + (11//7)w^2,
    ]
    scaled = scale_denominators(fs)
    @test scaled == [
                    2xz^2 + 6yz^2 + 5zz,
                    change_base_ring(
                        ZZ,
                        210*( (3//2)x^2 + (5//3)y^2 + (7//5)z^2 + (11//7)w^2 ),
                        parent=parent(first(scaled)))
                    ]


end

@testset "Modular reduction" begin
    modulo = 5
    fs = [
        1xz, 2yz, 8zz, 10wz,
        4xz + 5yz + 6yz^2
    ]
    reduced = reduce_modulo(fs, modulo)
    Rgf, (xg, yg, zg, wg) = parent(first(reduced)), gens(parent(first(reduced)))
    @test reduced == [xg, 2yg, 3zg, 0wg, 4xg + yg^2]

    f = 5xz
    freduced = reduce_modulo(f, modulo)
    @test freduced == [0]

end

@testset "Modular reduction types" begin


end


@testset "Rational reconstruction" begin
    ms = [107, 199, 509, 2^31 - 1]
    as = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    for m in ms
        for a in as
            gf = GF(m)
            ac = lift( numerator(a) // gf(denominator(a)) )
            ar1 = rational_reconstruction(ac, BigInt(m))
            ar2 = rational_reconstruction(Int(ac), m)
            @test ar1 == ar2 == a
        end
    end

    # random tests
    ms = Primes.primes(2^20, 2^21)[1:3]
    as = [QQ(rand(1:100), rand(1:100)) for _ in 1:10]
    for m in ms
        for a in as
            gf = GF(m)
            ac = lift( numerator(a) // gf(denominator(a)) )
            ar1 = rational_reconstruction(ac, BigInt(m))
            ar2 = rational_reconstruction(Int(ac), m)
            @test ar1 == ar2 == a
        end
    end

end



@testset "Chineese remainder theorem" begin

    modular_images(a, ms) = map(m -> mod(a, m), ms)

    ms = [2, 3, 5, 7]
    P = prod(ms)

    nums = [rand(1:P-1) for _ in 1:10]

    for a in nums
        rs = modular_images(a, ms)
        acrt = crt(rs, ms)
        @test acrt == a
    end

end
