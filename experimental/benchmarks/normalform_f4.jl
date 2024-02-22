
include("../../src/Groebner.jl")

if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using Nemo
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 4

system = Groebner.cyclicn(7, k=Nemo.QQ, internal_ordering=:degrevlex)
gb = Groebner.groebner(system);
xs = gens(parent(system[1]))

@time Groebner.normalform(gb, system, check=false);
collect([AbstractAlgebra.normal_form(s, gb) for s in system])

@benchmark Groebner.normalform($gb, $system, check=false)
@benchmark for s in system
    AbstractAlgebra.normal_form(s, gb)
end
