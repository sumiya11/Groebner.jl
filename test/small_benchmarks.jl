
@testset "qq small benchmarks" begin

    systems = [
        ("noon 5", Groebner.change_ordering(Groebner.noonn(5), :degrevlex)),
        ("noon 6", Groebner.change_ordering(Groebner.noonn(6), :degrevlex)),
        ("noon 7", Groebner.change_ordering(Groebner.noonn(7), :degrevlex)),
        ("cyclic 6", Groebner.change_ordering(Groebner.cyclicn(6), :degrevlex)),
        # ("cyclic 7", Groebner.change_ordering(Groebner.cyclicn(7), :degrevlex)),
        ("katsura 5", Groebner.change_ordering(Groebner.katsuran(5), :degrevlex)),
        ("katsura 6", Groebner.change_ordering(Groebner.katsuran(6), :degrevlex)),
        ("katsura 7", Groebner.change_ordering(Groebner.katsuran(7), :degrevlex)),
        ("henrion 5", Groebner.change_ordering(Groebner.henrion5(), :degrevlex)),
        # ("henrion 6", Groebner.change_ordering(Groebner.henrion6(), :degrevlex)),
    ]

    for (name, system) in systems
        @info "testing" name
        @test Groebner.isgroebner(Groebner.groebner(system))
    end

end

