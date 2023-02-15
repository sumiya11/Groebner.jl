
import Nemo
using Random

orderings_to_test = [Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
grounds_to_test = [Nemo.GF(2^31-1), Nemo.QQ]
representations_to_test = [
    Groebner.default_safe_representation(Groebner.NotPacked{UInt64}()),
    Groebner.Representation{Groebner.PowerVector{UInt64}}(),
    Groebner.Representation{Groebner.PowerVector{UInt32}}(),
    # Groebner.Representation{Groebner.PackedVector{UInt64, UInt8}}(),
    # Groebner.Representation{Groebner.PackedVector{UInt32, UInt16}}(),
    Groebner.Representation{Groebner.PackedPair1{UInt64, UInt8}}(),
    Groebner.Representation{Groebner.PackedPair2{UInt64, UInt8}}(),
    Groebner.Representation{Groebner.PackedPair3{UInt64, UInt8}}()
]

@testset "input-output nemo" begin
    rng = Random.MersenneTwister(42)

    for ord in orderings_to_test
        for ground in grounds_to_test
            for representation in representations_to_test
                ord_sym = Groebner.ordering_typed2sym(ord)
                R, (x, y) = Nemo.PolynomialRing(ground, ["x", "y"], ordering=ord_sym)
                fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
                ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
                meta = Groebner.set_metaparameters(ring, ord, false, false, :exact, rng)
                fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
                @test fsfs == fs

                root = Groebner.change_ordering(
                    Groebner.rootn(6, ground=ground, np=Nemo), 
                    ord_sym
                )
                ring, exps, cfs = Groebner.convert_to_internal(representation, root, Groebner.InputOrdering())
                meta = Groebner.set_metaparameters(ring, ord, false, false, :exact, rng)
                fsfs = Groebner.convert_to_output(ring, root, exps, cfs, meta)
                @test fsfs == root

                noon = Groebner.change_ordering(
                    Groebner.noonn(3, ground=ground, np=Nemo), 
                    ord_sym
                )
                ring, exps, cfs = Groebner.convert_to_internal(representation, noon, Groebner.InputOrdering())
                meta = Groebner.set_metaparameters(ring, ord, false, false, :exact, rng)
                fsfs = Groebner.convert_to_output(ring, noon, exps, cfs, meta)
                @test fsfs == noon

                for nn in (5, 10, 25)
                    R, xs = Nemo.PolynomialRing(ground, ["x$i" for i in 1:nn], ordering=ord_sym)
                    if Groebner.capacity(representation) >= nn
                        ring, exps, cfs = Groebner.convert_to_internal(representation, xs, Groebner.InputOrdering())
                        meta = Groebner.set_metaparameters(ring, ord, false, false, :exact, rng)
                        xsxs = Groebner.convert_to_output(ring, xs, exps, cfs, meta)
                        @test xsxs == xs
                    else
                        @test_throws AssertionError Groebner.convert_to_internal(representation, xs, Groebner.InputOrdering())
                    end
                end
            end
        end
    end
end

@testset "input-output consistency nemo" begin
    for ord in orderings_to_test
        for ground in grounds_to_test
            ord_sym = Groebner.ordering_typed2sym(ord)
            R, (x, y) = Nemo.PolynomialRing(ground, ["x", "y"], ordering=ord_sym)
            fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
        end
    end
end

@testset "input-output generic nemo" begin

    @test Groebner.determinechartype(2) === UInt8
    @test Groebner.determinechartype(2^5) === UInt16
    @test Groebner.determinechartype(2^10) === UInt32
    @test Groebner.determinechartype(2^31-1) === UInt64
    @test Groebner.determinechartype(2^60) === UInt128

    rng = Random.MersenneTwister(42)

    for ord in orderings_to_test
        for ground in grounds_to_test
            for representation in representations_to_test
                ord_sym = Groebner.ordering_typed2sym(ord)
                R, (x, y) = Nemo.PolynomialRing(ground, ["x", "y"], ordering=ord_sym)
                fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
                ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
                ring.origring = :hasparent
                meta = Groebner.set_metaparameters(ring, ord, false, false, :exact, rng)
                fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
                @test fsfs == fs
            end
        end
    end
end
