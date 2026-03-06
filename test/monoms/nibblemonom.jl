using Test, Groebner, Random

@testset "NibbleMonom, NibbleNoDeg" begin
    for NM in [Groebner.NibbleMonom{16}, Groebner.NibbleNoDeg{16}]

        x = [1, 2, 3, 0, 4]
        ev = Groebner.monom_construct_from_vector(NM, x)
        @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
        @test Groebner.monom_totaldeg(ev) === UInt64(10)
        @test Groebner.monom_entrytype(ev) === UInt8
        tmp = similar(x)
        @test Groebner.monom_to_vector!(tmp, ev) == [1, 2, 3, 0, 4]
        @test tmp == [1, 2, 3, 0, 4]

        @test Groebner.monom_max_vars(NM) == 32

        y = [10, 10, 10, 16, 10]
        @test_throws Groebner.MonomialDegreeOverflow Groebner.monom_construct_from_vector(
            NM,
            y
        )
    end

    N = 512
    NM = Groebner.NibbleNoDeg{N}
    x = [15 for i in 1:2N]
    m = Groebner.monom_construct_from_vector(NM, x)
    @test Groebner.monom_totaldeg(m) == sum(x)

    for i in 1:10
        N = 2^i
        for NM in [Groebner.NibbleMonom{N}, Groebner.NibbleNoDeg{N}]
        for k in 1:100
            for (x1,x2) in [
                (rand(0:15, 2N), rand(0:15, 2N)),
                (rand(0:1, 2N), rand(0:1, 2N)),
                (vcat([1,2,3], zeros(Int, 2N-3)), vcat([4,5,6], zeros(Int, 2N-3))),
            ]
                if NM == Groebner.NibbleMonom{N} && sum(x1) + sum(x2) > typemax(Groebner.monom_entrytype(NM)) 
                    continue
                end
                if any(>(15), x1 .+ x2) continue end
                m1 = Groebner.monom_construct_from_vector(NM, x1)
                m2 = Groebner.monom_construct_from_vector(NM, x2)
                m3 = Groebner.monom_product!(m1, m1, m2)
                x3 = Groebner.monom_to_vector!(similar(x1), m3)
                @test x3 == (x1 .+ x2)
                flag, m4 = Groebner.monom_is_divisible!(m1, m3, m2)
                @test flag && Groebner.monom_is_equal(m4, m1)
                m5 = Groebner.monom_lcm!(m4, m1, m2)
                x5 = Groebner.monom_to_vector!(similar(x1), m5)
                @test x5 == max.(x1, x2)
                flag = Groebner.monom_is_gcd_const(m1, m2)
                @test flag == all(iszero, x1 .* x2)
                rng = Random.Xoshiro(42)
                hv = Groebner.monom_construct_hash_vector(rng, NM, 2N)
                h1, h2 = Groebner.monom_hash(m1, hv), Groebner.monom_hash(m2, hv)
                @test h1 + h2 == Groebner.monom_hash(m3, hv)
            end
        end
        end
    end
end
