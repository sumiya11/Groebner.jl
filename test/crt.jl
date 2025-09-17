using Test, Groebner

@testset "crt" begin
    cases = []

    push!(
        cases,
        Dict(
            :moduli => Vector{UInt64}([1099511627791, 3518437208889]),
            :ai => Vector{UInt64}([8590035092, 8589936485])
        ),
        Dict(
            :moduli => Vector{UInt64}([1099511627791, 17, 3518437208889]),
            :ai => Vector{UInt64}([8590035092, 3, 8589936485])
        ),
        Dict(
            :moduli => Vector{UInt64}([1099511627791, 2^30 + 3, 3518437208889]),
            :ai => Vector{UInt64}([8590035092, 42, 0])
        ),
        Dict(
            :moduli => Vector{UInt64}([1099511627791, 2^30 + 3, 3518437208889]),
            :ai => Vector{UInt64}([0, 0, 0])
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
        @test 0 <= buf < M
        for i in 1:length(ci)
            @test (buf - ai[i]) % moduli[i] == 0
        end
    end
end

@testset "crt table" begin
    table_zz = [BigInt[], [BigInt(), BigInt()]]
    moduli = UInt64.([2, 3])
    res = [BigInt[], [BigInt(5), BigInt(0)]]
    modulo = BigInt()
    tables_ff = [[Int32[mod(res[i][j], moduli[k]) for j in 1:length(table_zz[i])] for i in 1:length(table_zz)] for k in 1:length(moduli)]
    mask = [falses(length(table_zz[i])) for i in 1:length(table_zz)]
    Groebner.crt_vec_full!(table_zz, modulo, tables_ff, moduli, mask)
    @test table_zz == res

    table_zz = [[BigInt(), BigInt(), BigInt()], [BigInt()], [BigInt(), BigInt()]]
    moduli = UInt64.([31, 39, 41])
    res = [[BigInt(prod(moduli) - 1), BigInt(21), BigInt(39)], [BigInt(1)], [BigInt(0), BigInt(7)]]
    modulo = BigInt()
    tables_ff = [[[mod(res[i][j], moduli[k]) for j in 1:length(table_zz[i])] for i in 1:length(table_zz)] for k in 1:length(moduli)]
    mask = [falses(length(table_zz[i])) for i in 1:length(table_zz)]
    for tasks in [1,2,3,4]
        Groebner.crt_vec_full!(table_zz, modulo, tables_ff, moduli, mask; tasks=tasks)
        @test table_zz == res
        @test modulo == prod(BigInt, moduli)
    end
end
