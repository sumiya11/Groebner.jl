# Use either of:
include((@__DIR__) * "/../../src/Groebner.jl")
# using Groebner

import AbstractAlgebra

using Base.Threads
import Pkg
using Logging
using BenchmarkTools

params = (linalg=:prob, threading=true)
ground = AbstractAlgebra.QQ
ord = :degrevlex

# compile
system = Groebner.katsuran(8, ground=ground, ordering=ord)
Groebner.groebner(system; loglevel=Logging.Error, params...)
system = Groebner.cyclicn(6, ground=ground, ordering=ord)
Groebner.groebner(system; loglevel=Logging.Error, params...)

@warn "Using $(nthreads()) threads."
@warn "Params: $params"

@info "cyclic 7.."
system = Groebner.cyclicn(7, ground=ground, ordering=ord)
@time Groebner.groebner(system; params...)

@info "katsura 8.."
system = Groebner.katsuran(8, ground=ground, ordering=ord)
@time Groebner.groebner(system; params...)

@info "katsura 9.."
system = Groebner.katsuran(9, ground=ground, ordering=ord)
@time Groebner.groebner(system; params...)

@info "eco 11.."
system = Groebner.eco11(ground=ground, ordering=ord)
@time Groebner.groebner(system; params...)

@info "usual corner case"
R, (x, y) = ground["x","y"]
system = [x^2*y + BigInt(10)^2000, x*y^2 + BigInt(11)^2000]
Groebner.groebner(system; params...)
@time Groebner.groebner(system; params...)
