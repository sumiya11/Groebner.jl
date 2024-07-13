import Primes

@testset "rational reconstruction" begin
    moduli =
        BigInt[107, 199, 509, 2^31 - 1, nextprime(BigInt(2)^100), nextprime(BigInt(2)^200)]
    numbers = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    for m in moduli
        bnd = Groebner.ratrec_reconstruction_bound(m)
        for a in numbers
            ac = mod(numerator(a) * invmod(denominator(a), m), m)
            success, pq = Groebner.ratrec_nemo(
                Groebner.Nemo.ZZRingElem(ac),
                Groebner.Nemo.ZZRingElem(m)
            )
            (num, den) = numerator(pq), denominator(pq)
            @test Base.unsafe_rational(num, den) == a
        end
    end

    # random checks
    moduli = BigInt[
        17,
        prod(Primes.nextprimes(BigInt(2^10), 20)),
        2^31 - 1,
        prod(Primes.nextprimes(BigInt(2^30 + 3), 5))
    ]
    samples = 1000
    for m in moduli
        bnd = Groebner.ratrec_reconstruction_bound(m)
        for i in 1:samples
            a = BigInt(rand(0:(m - 1)))
            ac = mod(numerator(a) * invmod(denominator(a), m), m)
            success, pq = Groebner.ratrec_nemo(
                Groebner.Nemo.ZZRingElem(ac),
                Groebner.Nemo.ZZRingElem(m)
            )
            (num, den) = numerator(pq), denominator(pq)
            if a < sqrt(div(m, 2))
                @test success
            end
            if success
                @test mod(num, m) == mod(a * den, m)
            end
        end
    end
end

function table_compute_res(table, m, T)
    [
        map(c -> T(mod(mod(numerator(c), m) * invmod(denominator(c), m), m)), table[j]) for
        j in 1:length(table)
    ]
end

@testset "vec ratrec + crt" begin
    table_qq = [[BigInt(0) // BigInt(1)], [BigInt(0) // BigInt(1), BigInt(0) // BigInt(1)]]
    table_zz = [[BigInt(58)], [BigInt(15), BigInt(73)]]
    modulo = BigInt(77)
    success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo)
    @assert success
    @assert table_qq[1][1] == BigInt(1) // BigInt(4)
    @assert table_qq[2][1] == BigInt(-2) // BigInt(5)
    @test table_qq[2][2] == BigInt(-4) // BigInt(1)

    # Impossible reconstruction
    table_qq = [[BigInt(0) // BigInt(1)]]
    table_zz = [[BigInt(643465418)]]
    modulo = BigInt(2^31 - 1)
    indices = [(1, 1)]
    success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo)
    @test !success

    boot = 10
    for _ in 1:boot
        k = rand(1:100)
        moduli = map(UInt32, Primes.nextprimes(rand((2^20):(2^28)), k))
        m = rand(1:10)
        P0 = prod(map(BigInt, moduli))
        P1 = floor(BigInt, sqrt(P0 / 2) - 1)
        N, D = P1, P1
        _s = rand(1:10, m)
        # Generate rational numbers
        table_ans = [rand((-N):N, _s[_m]) .// rand(1:D, _s[_m]) for _m in 1:m]
        # Compute k tables of residuals
        tables_res = [table_compute_res(table_ans, moduli[i], UInt32) for i in 1:k]
        # Run CRT
        table_zz =
            [[BigInt(0) for _ in 1:length(table_ans[i])] for i in 1:length(table_ans)]
        Groebner.crt_vec_full!(table_zz, modulo, tables_res, moduli)
        @test modulo == prod(BigInt, moduli)
        @test table_zz == table_compute_res(table_ans, modulo, BigInt)
        # Run rat. rec.
        table_qq = [
            [Rational{BigInt}(0) for _ in 1:length(table_ans[i])] for
            i in 1:length(table_ans)
        ]
        success = Groebner.ratrec_vec_full!(table_qq, table_zz, modulo)
        @test table_ans == table_qq && success
    end
end
