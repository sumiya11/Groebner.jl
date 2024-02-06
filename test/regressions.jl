import Primes

@testset "regression, SI.jl normalform" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])

    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)
    @test Groebner.normalform([x], [R(0)]) == [R(0)]
    @test Groebner.normalform([x], [R(1)]) == [R(1)]
    @test Groebner.normalform([x], [R(0), R(1), R(0)]) == [R(0), R(1), R(0)]
end

@testset "regression, ordering of empty" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], ordering=:lex)

    ord = Groebner.DegRevLex()
    gb1 = Groebner.groebner([x, y], ordering=ord)
    gb2 = Groebner.groebner([R(0)], ordering=ord)
    @test parent(first(gb1)) == R
    @test parent(first(gb2)) == R
    nf1 = Groebner.normalform([x, y], x, ordering=ord)
    nf2 = Groebner.normalform([x, y], R(0), ordering=ord)
    @test parent(nf1) == R
    @test parent(nf2) == R
end

@testset "regression, SI.jl cmp" begin
    # this may crash if the comparator is invalid
    function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}
        inds = collect(1:length(gens))
        cmp  = (x, y) -> Groebner.monom_isless(exps[gens[x]], exps[gens[y]], Groebner._DegRevLex{true}([]))
        sort!(inds, lt=cmp)
        @test true
    end

    gens = [2,2,2,2,2,2,2,2,2,
        3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        3,2,2,2,2,2,2,2,2,2,2,2,2,2,
        3,2,2,2,2,2
    ]

    exps = Groebner.ExponentVector{UInt64}[
        [0x0000, 0x0000],
        [0x0002, 0x0002],
        [0x0002, 0x0002]
    ]

    bug(gens, exps)
end

@testset "regression, column order in normalform" begin
    P, (x1, x2) = polynomial_ring(GF(1119649), ["x1", "x2"])

    gb = [x1 + 1]
    f = x1 + x2^2
    f_nf = Groebner.normalform(gb, f, ordering=Groebner.DegRevLex())
    residual = Groebner.normalform(gb, f - f_nf, ordering=Groebner.DegRevLex())

    @test iszero(residual)
end

@testset "regression, tracing invariants" begin
    R, (x1,x2,x3,x4,x5,x6,x7,_Z) = polynomial_ring(ZZ, [:x1,:x2,:x3,:x4,:x5,:x6,:x7,:_Z],ordering=:degrevlex)
    sys_z_t = [2*x1^2 - 2*x2^2 + 2*x3^2 - 2*x4^2 + 2*x5^2 - 2*x6^2 + 2*x7^2 - 1
        2*x1^3 - 2*x2^3 + 2*x3^3 - 2*x4^3 + 2*x5^3 - 2*x6^3 + 2*x7^3 - 1
        2*x1^4 - 2*x2^4 + 2*x3^4 - 2*x4^4 + 2*x5^4 - 2*x6^4 + 2*x7^4 - 1
        2*x1^5 - 2*x2^5 + 2*x3^5 - 2*x4^5 + 2*x5^5 - 2*x6^5 + 2*x7^5 - 1
        2*x1^6 - 2*x2^6 + 2*x3^6 - 2*x4^6 + 2*x5^6 - 2*x6^6 + 2*x7^6 - 1
        2*x1^7 - 2*x2^7 + 2*x3^7 - 2*x4^7 + 2*x5^7 - 2*x6^7 + 2*x7^7 - 1
        2*x1^8 - 2*x2^8 + 2*x3^8 - 2*x4^8 + 2*x5^8 - 2*x6^8 + 2*x7^8 - 1
        -35*x1 + 43*x2 + 63*x3 - 80*x4 - 56*x5 + 73*x6 - 73*x7 + _Z]
    pr = 2^28
    bad = 268419493
    sys_zp0 = map(f -> AbstractAlgebra.change_base_ring(AbstractAlgebra.GF(Primes.prevprime(pr-1)), f), sys_z_t)
    sys_zp2 = map(f -> AbstractAlgebra.change_base_ring(AbstractAlgebra.GF(bad), f), sys_z_t)
    trace, _ = Groebner.groebner_learn(sys_zp0, ordering=Groebner.DegRevLex());
    flag,gb = Groebner.groebner_apply!(trace, sys_zp2);
    @test true
end
