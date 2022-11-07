
function n_variable_set(n)
    R, x = PolynomialRing(QQ, ["x$i" for i in 1:n]);
    f = [
        sum(prod(x[i:n-k], init=1) for i in 1:k+1)
        for k in 0:n-1
    ];
    f
end

function test_n_variables(n)
    f = n_variable_set(n)
    x = gens(parent(f[1]))
    gb = Groebner.groebner(f)

    evencf(i) = isone(i) ? 0//1 : (2(i-1))//(2i-1)
    oddcf(i) = isone(i) ? 0//1 : (2(i-1)-1)//(2(i-1))
    cf(i, n) = (iseven(n) ? oddcf(i) : evencf(i))*(-1n)^(i == div(n, 2)+1)

    ans = [
        x[div(n, 2) - i + 2] - cf(i, n)
        for i in 1:div(n, 2)+1
    ]
    
    @test Groebner.isgroebner(gb)
    @test all(iszero, Groebner.normalform(gb, f))
    @test gb == ans
end

@testset "handling many variables" begin
    
    # up to 63
    for n in 2:3:63
        test_n_variables(n)
    end

    # up to 127
    for n in [64, 100, 101, 127]
        @info "" n
        test_n_variables(n)
    end

    # up to 511
    for n in [128, 257, 256, 511]
        @info "" n
        R, x = PolynomialRing(QQ, ["x$i" for i in 1:n])
        f = x
        Groebner.groebner(f)
    end

end
