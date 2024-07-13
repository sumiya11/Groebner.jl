using AbstractAlgebra
using Primes

function get_test_system1(R, N)
    (x1,x2,x3,x4) = AbstractAlgebra.gens(R)
    system = [
        x1 + x2 + x3 + x4,
        x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
        x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
        x1 * x2 * x3 * x4 + N
    ]
    result = [
        x1 + x2 + x3 + x4,
        x2^2 + x2 * x3 + x3^2 + x2 * x4 + x3 * x4 + x4^2,
        x3^3 + x3^2 * x4 + x3 * x4^2 + x4^3,
        x4^4 - N
    ]
    system, result
end

@testset "groebner modular-hard problems" begin
    R, (x1, x2, x3, x4) =
        polynomial_ring(QQ, ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)

    N = prod(map(BigInt, nextprimes(2^30 + 3, 5)))
    system, result = get_test_system1(R, N)
    # this should take about 10 primes
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 5)))
    system, result = get_test_system1(R, N)
    # this should take about 10 primes
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 100)))
    system, result = get_test_system1(R, N)
    # around 200 primes are required
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 5_000)))
    # around 10k primes are required
    system, result = get_test_system1(R, N)
    gb = Groebner.groebner(system)
    @test gb == result

    # TODO: Sasha is too greedy to support this case.
    # system = [x1 - (2^31 - 1) * x2 - (2^30 + 3) * x3]
    # @test Groebner.groebner(system) == system

    for start_of_range in [2^10, 2^20, 2^30, 2^40]
        for size_of_range in [10, 100, 1000]
            N = prod(Primes.nextprimes(BigInt(start_of_range), size_of_range))
            system = [x1 + N // (N + 1), x3 - (N + 1) // (N - 1), x2 + 42]
            gb = Groebner.groebner(system)
            @test gb == [x3 - (N + 1) // (N - 1), x2 + 42, x1 + N // (N + 1)]
        end
    end
end

@testset "groebner strange poly" begin
    get_rand_poly(x, d, n) = sum([rand(1:5) * prod(x .^ rand(0:d, length(x))) for _ in 1:n])
    
    n = 10
    R, x = polynomial_ring(GF(17), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    f1 = get_rand_poly(x, 10, 2^15);
    @test Groebner.groebner([f1]) == [divexact(f1, leading_coefficient(f1))]

    n = 6
    R, x = polynomial_ring(GF(17), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    f1 = get_rand_poly(x, 3, 2^10)
    f2 = get_rand_poly(x, 2, 2^10)
    gb = Groebner.groebner([f1, f2])
    @info "GB contains polynomials of lengths: $(sort(map(length, gb)))"
end
