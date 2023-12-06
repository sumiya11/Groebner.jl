using AbstractAlgebra, AllocCheck, BenchmarkTools, Logging, IOCapture

function get_all_monoms_up_to_total_degree(ring, vars, degree)
    if degree == 0 || length(vars) == 0
        return [one(ring)]
    end
    res = get_all_monoms_up_to_total_degree(ring, vars, degree - 1)
    res_prev = copy(res)
    for var in vars
        append!(res, res_prev .* var)
    end
    unique!(res)
end

R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)
s = [2x * y + 5y + 7, 9x * y + 11x + 13]

g = GF(2^31 - 1)
o = :degrevlex
si = Groebner.eco5(ordering=:degrevlex, ground=GF(2^31 - 1));

# R, x = PolynomialRing(GF(2^31 - 1), [["x$i" for i in 1:6]...], ordering=:degrevlex)
# m = get_all_monoms_up_to_total_degree(R, gens(R), 5);
# length(m)
# S = map(i -> sum(rand(m, div(length(m), 2))), 1:7);
# map(length, S)
# @time Groebner.groebner(S, linalg=:deterministic);
# si = S
graph, gb = Groebner.groebner_learn(si);
@time flag, gb2 = Groebner.groebner_apply!(graph, si, loglevel=0);
@assert flag

begin
    io1 = open((@__DIR__) * "/logs_direct.txt", "w")
    io2 = open((@__DIR__) * "/logs_deterministic.txt", "w")
    io3 = open((@__DIR__) * "/logs_randomized.txt", "w")
    io4 = open((@__DIR__) * "/logs_apply.txt", "w")

    c1 = IOCapture.capture() do
        @time gb2 =
            Groebner.groebner(si, linalg=:direct_rref, sparsity=:sparse, loglevel=-3)
    end
    println(io1, c1.output)

    c2 = IOCapture.capture() do
        @time gb2 = Groebner.groebner(si, linalg=:deterministic, loglevel=-3)
    end
    println(io2, c2.output)

    c3 = IOCapture.capture() do
        @time gb2 = Groebner.groebner(si, linalg=:randomized, loglevel=-3)
    end
    println(io3, c3.output)

    graph, gb = Groebner.groebner_learn(si)
    c4 = IOCapture.capture() do
        @time flag, gb2 = Groebner.groebner_apply!(graph, si, loglevel=-3)
        @assert flag
        gb2
    end
    println(io4, c4.output)

    @assert c1.value == c2.value == c3.value == c4.value

    close(io1)
    close(io2)
    close(io3)
    close(io4)
end

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
