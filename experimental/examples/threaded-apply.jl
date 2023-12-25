using Groebner, AbstractAlgebra, Primes, Base.Threads
using BenchmarkTools

Groebner.logging_enabled() = false

# Computes the bases of the given system modulo different primes
function compute_bases(system, iters::Int, batch_size::Int)
    prime = 2^30 + 3
    Zp = GF(prime)

    system_zp = map(f -> map_coefficients(c -> Zp(c), f), system)
    trace, _ = groebner_learn(system_zp)

    bases = Vector{typeof(system_zp)}()
    for i in 1:iters
        for j in 1:batch_size
            prime = nextprime(prime + 1)
            Zp = GF(prime)
            system_zp = map(f -> map_coefficients(c -> Zp(c), f), system)

            flag, gb = groebner_apply!(trace, system_zp)
            @assert flag

            push!(bases, gb)
        end
    end

    return bases
end

# Same as above, but uses multi-threading
function compute_bases_threaded(system, iters::Int, batch_size::Int)
    prime = 2^30 + 3
    Zp = GF(prime)

    system_zp = map(f -> map_coefficients(c -> Zp(c), f), system)
    trace, _ = groebner_learn(system_zp)

    bases = Vector{typeof(system_zp)}()
    buffer_traces = map(_ -> deepcopy(trace), 1:nthreads())
    buffer_bases = map(_ -> empty(system_zp), 1:batch_size)

    for i in 1:iters
        buffer_primes = nextprimes(prime + 1, batch_size)

        # :static guarantees that threadid() persists within one iteraton
        Base.Threads.@threads :static for j in 1:batch_size
            prime = buffer_primes[j]
            Zp = GF(prime)
            system_zp = map(f -> map_coefficients(c -> Zp(c), f), system)

            threadlocal_trace = buffer_traces[threadid()]

            flag, gb = groebner_apply!(threadlocal_trace, system_zp)
            @assert flag

            buffer_bases[j] = gb
        end

        prime = last(buffer_primes)
        append!(bases, buffer_bases)
    end

    return bases
end

system = Groebner.katsuran(9, ordering=:degrevlex, ground=ZZ)
bases_1 = compute_bases(system, 5, 8);
bases_2 = compute_bases_threaded(system, 5, 8);
@assert bases_1 == bases_2

@btime compute_bases($system, 5, 8);
@btime compute_bases_threaded($system, 5, 8);
