using Random, BenchmarkTools

include("common.jl")
include("hashtable-1.jl")
include("hashtable-2.jl")

n = 10
sz = 2^16

function setup_random_monoms(n, d, s; T=UInt8)
    monoms = Vector{Vector{T}}(undef, s)
    for i in 1:s
        monoms[i] = rand(T(0):T(d), n)
    end
    monoms
end

function setup_1(n, sz, s)
    ht = hashtable_initialize1(n, Random.MersenneTwister(42), Vector{UInt8}, sz)
    monoms = setup_random_monoms(n, 5, s)
    ht, monoms
end

function setup_2(n, sz, s)
    ht = hashtable_initialize1(n, Random.MersenneTwister(42), Vector{UInt8}, sz)
    monoms = setup_random_monoms(n, 5, s)
    for i in 1:length(monoms)
        hashtable_insert!(ht, monoms[i])
    end
    monoms2 = setup_random_monoms(n, 5, s)
    for j in 1:length(monoms)
        if iszero(j % 1_000)
            monoms[j] = monoms2[j]
        end
    end
    ht, monoms
end

# Inserts, almost no collisions
for k in 8:16
    sz = 2^k
    s = div(sz, 2)
    @info "n = $n, sz = 2^$k, s = $s"

    @btime begin
        for i in 1:($s)
            hashtable_insert!(ht, monoms[i])
        end
    end setup = begin
        ht, monoms = setup_1($n, $sz, $s)
    end
end

# Inserts, almost all are hits
for k in 8:16
    sz = 2^k
    s = div(sz, 2)
    @info "n = $n, sz = 2^$k, s = $s"

    @btime begin
        for i in 1:($s)
            hashtable_insert!(ht, monoms[i])
        end
    end setup = begin
        ht, monoms = setup_2($n, $sz, $s)
    end
end
