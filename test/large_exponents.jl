
# We guarantee correctness up to total degrees of 2^32-1

@testset "large exponents handling" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex) 

    # up to 2^8-1
    for (i, d) in enumerate(4:2:255)
        # i = [...] are not in the sequence
        d in [124, 30] && continue
        f = [x^d - 1, x*y + 2]
        m, n, k = div(d, 2), div(d, 2) + 1, div(d, 2) - 1
        gb = Groebner.groebner(f)
        @test gb == [
            x*y + 2, 
            x^m + (-1)^(i)//BigInt(2)^m*y^m, 
            y^n + (-1)^(i+1)*BigInt(2)^n*x^k
        ]
    end

    # up to 5^6 < 2^14
    for i in 1:6
        u, v = 3^i, 5^i
        f = [x^u*y^v - 1, x^v + y^u]
        gb = Groebner.groebner(f)
        @test gb == [
            x^v + y^u,
            y^(v + u) + x^(v - u),
            x^u*y^v - 1
        ]
    end

    # above 2^16-1
    f = [x^(2^16) + y]
    @test f == Groebner.groebner(f)

    # up to 5^13 < 2^32
    for i in 7:13
        u, v = 3^i, 5^i
        f = [x^u*y^v - 1, x^v + y^u]
        gb = Groebner.groebner(f)
        @test gb == [
            x^v + y^u,
            y^(v + u) + x^(v - u),
            x^u*y^v - 1
        ]
    end

    # 2^30
    f = [x^1073741824 + y]
    @test f == Groebner.groebner(f)

end
