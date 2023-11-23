using AbstractAlgebra, AllocCheck, BenchmarkTools

R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)

s = [2x * y + 5y + 7, 9x * y + 11x + 13]

g = GF(2^31 - 1)
o = :degrevlex
si = Groebner.katsuran(8, ordering=:degrevlex, ground=GF(2^31 - 1))

@time gb1 = Groebner.groebner(si, linalg=:direct_rref, sparsity=:sparsedense, loglevel=0);

@time gb2 = Groebner.groebner(si, linalg=:deterministic);

for si in [
    Groebner.katsuran(6, ordering=o, ground=g),
    Groebner.katsuran(7, ordering=o, ground=g),
    Groebner.katsuran(8, ordering=o, ground=g),
    Groebner.katsuran(9, ordering=o, ground=g),
    Groebner.noonn(7, ordering=o, ground=g),
    Groebner.cyclicn(6, ordering=o, ground=g),
    Groebner.cyclicn(7, ordering=o, ground=g),
    Groebner.rootn(7, ordering=o, ground=g),
    Groebner.rootn(8, ordering=o, ground=g),
    Groebner.reimern(6, ordering=o, ground=g),
    Groebner.eco7(ordering=o, ground=g)
]
    println("###")
    GC.gc()
    @time gb1 =
        Groebner.groebner(si, linalg=:direct_rref, sparsity=:sparsedense, loglevel=0)
    println(Groebner.count1[] / Groebner.count2[])
    @time gb2 = Groebner.groebner(si, linalg=:direct_rref, sparsity=:sparse, loglevel=0)
    @time gb3 = Groebner.groebner(si, linalg=:deterministic)
    @assert gb1 == gb2 == gb3
end

@btime Groebner.groebner($si, linalg=:direct_rref, sparsity=:sparsedense);
@btime Groebner.groebner($si, linalg=:deterministic);

@profview Groebner.groebner(si, linalg=:direct_rref, sparsity=:sparsedense);

gb1 == gb2
