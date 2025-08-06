using Groebner, Random, AbstractAlgebra, Primes

systems = [
    ("chandra-9", Groebner.Examples.chandran(9)),
    ("chandra-10", Groebner.Examples.chandran(10)),
    ("cyclic-7", Groebner.Examples.cyclicn(7)),
    ("eco-11", Groebner.Examples.eco11()),
    ("katsura-9", Groebner.Examples.katsuran(9)),
    ("hexapod", Groebner.Examples.hexapod()),
    ("hiv2", Groebner.Examples.HIV2()),
    ("alea6", Groebner.Examples.alea6())
]

if false
    for (name, sys) in systems
        @info "Running $name.."
        gb = groebner(sys)
        @time res1 = isgroebner(gb, linalg=:deterministic)
        @time res2 = isgroebner(gb, linalg=:randomized)
        @assert res1 == res2
    end
end

rng = Random.Xoshiro(42)

for k in [GF(prevprime(600)), GF(prevprime(2^10)), GF(prevprime(2^20))]
    boot = 100000
    for _ in 1:boot
        sys = Groebner.Examples.random_generating_set(
            rng,
            k,
            rand(rng, (:lex, :deglex, :degrevlex)),
            rand(rng, 2:4),
            rand(rng, 1:10),
            rand(rng, 1:10),
            rand(rng, 1:10),
            rand(rng, 100:(2 ^ 30))
        )
        if isempty(sys)
            continue
        end
        res1 = isgroebner(sys, linalg=:deterministic)
        res2 = isgroebner(sys, linalg=:randomized)
        if res1 != res2
            @info "Hehe"
        end
    end
end
