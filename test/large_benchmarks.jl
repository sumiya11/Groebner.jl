
ground = GF(2^31 - 1)
systems = [
    ("henrion 7", Groebner.henrion7(ground=ground)),
    ("noon 8", Groebner.noonn(8, ground=ground)),
    ("noon 9", Groebner.noonn(9, ground=ground)),
    ("cyclic 8", Groebner.cyclicn(8, ground=ground)),
    # ("cyclic 8", Groebner.cyclicn(9, ground=ground)),
    ("katsura 11", Groebner.katsuran(11, ground=ground)),
    # ("katsura 12", Groebner.katsuran(12)),
    # ("henrion 6", Groebner.henrion6()),
    ("eco 12", Groebner.eco12())
]

@testset "ff large benchmarks" begin
    for (name, system) in systems
        @info "testing" name
        system = Groebner.change_ordering(system, :degrevlex)
        @test Groebner.isgroebner(Groebner.groebner(system, linalg=:prob))
    end
end
