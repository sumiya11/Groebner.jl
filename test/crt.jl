@testset "CRT" begin
    cases = []

    push!(
        cases,
        Dict(
            :moduli => Vector{UInt64}([1099511627791, 3518437208889]),
            :ai => Vector{UInt64}([8590035092, 8589936485])
        )
    )

    for c in cases
        ai = c[:ai]
        moduli = c[:moduli]
        M = prod(Vector{BigInt}(moduli))
        ci = Vector{BigInt}([0 for _ in moduli])
        n1, n2 = BigInt(0), BigInt(0)
        Groebner.crt_precompute!(M, n1, n2, ci, moduli)
        for i in 1:length(ci)
            pii = div(M, moduli[i])
            @test ci[i] % pii == 0
            @test (ci[i] - 1) % moduli[i] == 0
        end
        buf = BigInt(0)
        Groebner.crt!(M, buf, n1, n2, ai, ci)
        @test buf < M
        for i in 1:length(ci)
            @test (buf - ai[i]) % moduli[i] == 0
        end
    end
end
