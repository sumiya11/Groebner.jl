using Revise, Groebner, Random, BenchmarkTools, Printf, SmallCollections

function construct_random_monom(MonomT, nvars, rng)
    ev = rand(rng, 0:3, nvars)
    return Groebner.monom_construct_from_vector(MonomT, ev)
end

function bench_insert_only(nvars, MonomT, n_insertions=500)
    ord = Groebner.DegRevLex()
    ring = Groebner.PolyRing(nvars, ord, Groebner.CoeffModular(1009), :zp)
    rng = Random.default_rng()
    
    monoms = [construct_random_monom(MonomT, nvars, rng) for _ in 1:n_insertions]
    
    function fill_many(ht, monoms)
        for m in monoms
            Groebner.hashtable_insert!(ht, m)
        end
    end

    # Case 1: Inserting into an empty hashtable
    b_not_filled = @benchmarkable $fill_many(ht, $monoms) setup=begin
        ht = Groebner.hashtable_initialize($ring, $rng, $MonomT);
        Groebner.hashtable_resize_if_needed!(ht, $(3*length(monoms)))
    end evals=1
    
    # Case 2: Inserting into an already filled hashtable
    ht_filled = Groebner.hashtable_initialize(ring, rng, MonomT)
    Groebner.hashtable_resize_if_needed!(ht_filled, 3*length(monoms))
    fill_many(ht_filled, monoms)
    
    b_filled = @benchmarkable $fill_many($ht_filled, $monoms)
    
    return b_not_filled, b_filled
end

function run_benchmarks(N)
    if N == 11
        configs = [
            (N, Groebner.PackedTuple2{UInt64, UInt8}),
            (N, Vector{UInt8}),
            (N, SmallCollections.FixedVector{nextpow(2,N), UInt8}),
            (N, Groebner.FixedMonom{nextpow(2,N), UInt8}),
            (N, Groebner.FixedMonomNoDeg{nextpow(2,N), UInt8}),
            (N, Groebner.NibbleMonom{nextpow(2,N) ÷ 2})
        ]
    elseif N == 55
        configs = [
            (N, Vector{UInt8}),
            (N, SmallCollections.FixedVector{nextpow(2,N), UInt8}),
            (N, Groebner.FixedMonom{nextpow(2,N), UInt8}),
            (N, Groebner.FixedMonomNoDeg{nextpow(2,N), UInt8}),
            (N, Groebner.NibbleMonom{nextpow(2,N) ÷ 2})
        ]
    end

    results = []
    
    for (nvars, MonomT) in configs
        insertions = 10000
        @info "Benchmarking $insertions insertions for $MonomT with nvars=$nvars"
        
        b_not_filled, b_filled = bench_insert_only(nvars, MonomT, insertions)
        
        res_not_filled = run(b_not_filled)
        res_filled = run(b_filled)
        
        t_not_filled = minimum(res_not_filled).time / insertions
        t_filled = minimum(res_filled).time / insertions
        
        monom_name = string(MonomT)
        monom_name = replace(monom_name, "Groebner." => "")
        
        push!(results, (monom_name, t_not_filled, t_filled))
    end
    
    println("\n=======================================================")
    println(" Benchmark Results for N = $N (Min ns / Insertion) ")
    println("=======================================================")
    @printf("%-30s | %-10s | %-10s\n", "Monomial Type", "Empty HT", "Filled HT")
    println("-------------------------------+------------+----------")
    for (name, t1, t2) in results
        @printf("%-30s | %8.2f   | %8.2f \n", name, t1, t2)
    end
    println("=======================================================\n")
end

run_benchmarks(11)
run_benchmarks(55)
