using BenchmarkTools
using Random
using Logging

using Groebner, AbstractAlgebra

const SUITE = BenchmarkGroup()

systems = Dict(
    "katsura-8" => Groebner.katsuran(8, ground=GF(2^31-1), ordering=:degrevlex),
    "cyclic-7" => Groebner.cyclicn(7, ground=GF(2^31-1), ordering=:degrevlex)
)

SUITE["groebner"] = BenchmarkGroup()
SUITE["groebner"]["finite-field"] = BenchmarkGroup()

SUITE["groebner"]["finite-field"]["katsura-8"] = @benchmarkable Groebner.groebner($systems["katsura-8"], linalg=:prob, loglevel=Logging.Warn)
SUITE["groebner"]["finite-field"]["cyclic-7"] = @benchmarkable Groebner.groebner($systems["cyclic-7"], linalg=:prob, loglevel=Logging.Warn)

