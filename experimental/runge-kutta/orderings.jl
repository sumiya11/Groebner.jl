
using AbstractAlgebra
import Combinatorics
using Plots
using Random
include("yang.jl")

R = parent(polys[1])

xs = gens(R)
times = Dict()
Groebner._perm[] = collect(1:length(xs))

@profview Groebner.groebner(polys, ordering=Groebner.DegRevLex(), linalg=:prob);

times[Groebner._perm[]] = t.time

i = 0
for ord in Combinatorics.permutations(collect(1:length(xs)))
    Groebner._perm[] = shuffle(collect(1:length(xs)))

    t = @timed Groebner.groebner(polys, ordering=Groebner.DegRevLex(), linalg=:prob)
    times[Groebner._perm[]] = t.time
    if i == 1
        break
    end
    i += 1
end
