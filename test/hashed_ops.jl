
using .GroebnerBases: HashedMonomial, HashedTerm, HashedPolynomial,
                convert_to_hash_repr, convert_to_original_repr,
                MonomialHashtable, add_monomial!,
                is_monom_divisible, monom_divide, term_divide,
                is_monom_gcd_constant, monom_gcd, monom_lcm,
                term_poly_product, monom_poly_product

@testset "HashedMonomial & HashedTerm" begin
    m1, m2 = HashedMonomial(8), HashedMonomial(UInt(5))
    @test hash(m1) == 8
    @test hash(m2) == 5
    @test m1 == HashedMonomial(8)

    t1 = HashedTerm{typeof(QQ(1))}(8, QQ(1))
    t2 = HashedTerm{typeof(GF(3)(3))}(0, GF(3)(2))
end

@testset "MPoly{GF} ~ hash repr" begin
    ground = GF(2^31 - 1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"])

    HT = MonomialHashtable(R)

    polys = [R(1), x, 4x^2, 3x*z + z*y, z^2*y - z*y, 2x*y*z, x*y^2*z + 3x, 3x^2]
    true_lengths = [1, 2, 3, 5, 6, 7, 8, 8]

    for (i, poly) in enumerate(polys)
        hashpoly = convert_to_hash_repr(poly, HT)
        @test hcat(map(m -> HT[m], hashpoly.monoms)...) == poly.exps
        @test hashpoly.coeffs == poly.coeffs
        @test hashpoly.parent == poly.parent

        origpoly = convert_to_original_repr(hashpoly, HT)
        @test origpoly == poly

        @test length(HT.expmap) == true_lengths[i]
    end
end

@testset "Hash monomial division" begin

    ground = GF(2^31-1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"])
    HT = MonomialHashtable(R)

    ms = [R(1), x, x, x^2, x*y, x*z, y^2, y*z, z^2,
            x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, x*z, z^2*y, x*y*z, x*y^2*z]

    for m1 in ms
        for m2 in Hs
            true_flag, true_q = divides(m1, m2)

            # test custom division components

            # is_term_divisible
            mh1 = leading_monomial(convert_to_hash_repr(m1, HT))
            mh2 = leading_monomial(convert_to_hash_repr(m2, HT))
            flag = is_monom_divisible(mh1, mh2, HT)
            @test flag == true_flag

            if flag
                q = monom_divide(mh1, mh2, HT)
                qpoly = HashedPolynomial([q], [ground(1)], R)
                @test convert_to_original_repr(qpoly, HT) == true_q
            end
        end
    end

end

@testset "Hash monom gcd and lcm" begin

    ground = GF(2^31-1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"])
    HT = MonomialHashtable(R)

    ms = [R(1), x, x, x^2, x*y, x*z, y^2, y*z, z^2,
            x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, x*z, z^2*y, x*y*z, x*y^2*z]

    for m1 in ms
        for m2 in Hs
            true_gcd = gcd(m1, m2)
            true_flag = isconstant(true_gcd)

            # test custom gcd components

            mh1 = leading_monomial(convert_to_hash_repr(m1, HT))
            mh2 = leading_monomial(convert_to_hash_repr(m2, HT))
            flag = is_monom_gcd_constant(mh1, mh2, HT)
            @test flag == true_flag

            q = monom_gcd(mh1, mh2, HT)
            qpoly = HashedPolynomial([q], [ground(1)], R)
            @test convert_to_original_repr(qpoly, HT) == true_gcd

            q = monom_lcm(mh1, mh2, HT)
            qpoly = HashedPolynomial([q], [ground(1)], R)
            @test convert_to_original_repr(qpoly, HT) == lcm(m1, m2)
        end
    end

end

@testset "Hash multiplication by monomial" begin

    ground = GF(2^31-1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"])
    HT = MonomialHashtable(R)

    ms = [R(1), x, 5x, x^2, x*y, 4x*z, y^2, y*z, z^2,
            3x*y*z, x^2*y^3*z^4]
    Hs = [R(1) + x, x, x^2 + z^2*y^2, x*z + 5, z^2*y - z*y, x*y*z,
          4x*y^2*z + x, 3x*y*z, x*y*z, 3x*y - 4]

    for m in ms
        for H in Hs
            true_prod = leading_monomial(m) * H
            mh1 = leading_monomial(convert_to_hash_repr(m, HT))
            Hh2 = convert_to_hash_repr(H, HT)

            my_prod = monom_poly_product(mh1, Hh2, HT)
            @test true_prod == convert_to_original_repr(my_prod, HT)

            ########

            true_prod = m * H

            mh1 = leading_term(convert_to_hash_repr(m, HT))
            my_prod = term_poly_product(mh1, Hh2, HT)
            @test true_prod == convert_to_original_repr(my_prod, HT)
        end
    end

end
