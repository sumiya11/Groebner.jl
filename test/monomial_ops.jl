
using .GroebnerBases: exponent_ordering, is_term_divisible,
                unsafe_term_divide, term_divides,
                is_term_gcd_constant, term_gcd, term_lcm,
                term_equal

@testset "AA helper functions" begin

    R, _ = PolynomialRing(QQ, ["x", "y"])
    @test exponent_ordering(R, 2) == 2:-1:1

    R, _ = PolynomialRing(QQ, ["x", "y", "z"], ordering=:degrevlex)
    @test exponent_ordering(R, 3) == 1:2

end

@testset "Custom term division" begin

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    ms = [R(1), x, 5x, x^2, x*y, 4x*z, (1//3)y^2, y*z, z^2,
            3x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, (1//8)x*z + 5, z^2*y - z*y, x*y*z, 4x*y^2*z + x]

    for m in ms
        for H in Hs
            true_flag, true_q = divides(m, leading_term(H))

            # test custom division components

            # is_term_divisible
            flag = is_term_divisible(m, H)
            @test flag == true_flag

            if flag
                q = unsafe_term_divide(m, H)
            else
                q = R(0)
            end
            @test q == true_q

            @test GroebnerBases.term_divides(m, H) == divides(m, leading_term(H))
        end
    end

end

@testset "Custom term gcd" begin

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    ms = [R(1), x, 5x, x^2, x*y, 4x*z, (1//3)y^2, y*z, z^2,
            3x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, (1//8)x*z + 5, z^2*y - z*y, x*y*z, 4x*y^2*z + x]

    for m in ms
        for H in Hs
            true_gcd = gcd(m, leading_term(H))
            my_gcd = GroebnerBases.term_gcd(m, H)

            @test true_gcd == my_gcd
            @test is_term_gcd_constant(m, H) == isconstant(gcd(m, leading_term(H)))
        end
    end

end


@testset "Custom term lcm" begin

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    ms = [R(1), x, 5x, x^2, x*y, 4x*z, (1//3)y^2, y*z, z^2,
            3x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, (1//8)x*z + 5, z^2*y - z*y, x*y*z, 4x*y^2*z + x]

    for m in ms
        for H in Hs
            true_lcm = lcm(m, leading_term(H))
            my_lcm = GroebnerBases.term_lcm(m, H)

            @test monomial(true_lcm, 1) == monomial(my_lcm, 1)
        end
    end

end

@testset "Custom term equal" begin

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    ms = [R(1), x, 5x, x^2, x*y, 4x*z, (1//3)y^2, y*z, z^2,
            3x*y*z, x^2*y^3*z^4]
    Hs = [R(1), x, x^2, (1//8)x*z + 5, z^2*y - z*y, x*y*z, 4x*y^2*z + x,
            3x*y*z, x*y*z]

    for m in ms
        for H in Hs
            true_eq = m == leading_term(H)
            my_eq = GroebnerBases.term_equal(m, H)

            @test true_eq == my_eq
        end
    end

end
