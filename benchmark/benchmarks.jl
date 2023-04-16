using BenchmarkTools
using Random
using Logging

using Groebner, AbstractAlgebra

const SUITE = BenchmarkGroup()

K = GF(2^31-1)
systems_over_ff = Dict(
    "katsura-8" => Groebner.katsuran(8, ground=K, ordering=:degrevlex),
    "katsura-9" => Groebner.katsuran(9, ground=K, ordering=:degrevlex),
    "cyclic-7" => Groebner.cyclicn(7, ground=K, ordering=:degrevlex),
    "cyclic-8" => Groebner.cyclicn(8, ground=K, ordering=:degrevlex),
    "noon-7"  => Groebner.noonn(7, ground=K, ordering=:degrevlex),
    "noon-8"  => Groebner.noonn(8, ground=K, ordering=:degrevlex),
    "eco-10"  => Groebner.eco10(ground=K, ordering=:degrevlex),
    "eco-11"  => Groebner.eco11(ground=K, ordering=:degrevlex)
)

run_system(system; kwargs...) = Groebner.groebner(system, linalg=:prob, kwargs...)

SUITE["groebner"] = BenchmarkGroup()
SUITE["groebner"]["finite-field"] = BenchmarkGroup()

for (name, system) in systems_over_ff
    SUITE["groebner"]["finite-field"][name] = @benchmarkable run_system($system)
end