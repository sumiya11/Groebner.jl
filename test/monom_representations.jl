
all_hint_types = [
    Groebner.NotPacked,
    Groebner.Packed
]

all_hint_params = [
    UInt8,
    UInt16,
    UInt32,
    UInt64
]

@testset "monom representation hints" begin
    T = UInt64
    repr = Groebner.Representation
    saferepr = Groebner.SafeRepresentation
    unsaferepr = Groebner.UnsafeRepresentation

    using AbstractAlgebra
    R, (x, y) = PolynomialRing(QQ, ["x","y"])
    f = [x + y]

    @test_logs (:warn,) Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.NotPacked{UInt8}())
    @test_logs (:warn,) Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.NotPacked{UInt16}())
    @test_logs Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.NotPacked{UInt32}())
    @test_logs Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.NotPacked{UInt64}())

    @test_logs (:warn,) Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.Packed{UInt8}())
    @test_logs (:warn,) Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.Packed{UInt16}())
    @test_logs Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.Packed{UInt32}())
    @test_logs (:warn,) Groebner.guess_effective_representation(f, saferepr(), Groebner.Lex(), Groebner.Packed{UInt64}())
    
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt8}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt16}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt32}())
    @test_logs (:warn,) Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt64}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt8}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt16}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt32}())
    @test_logs Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt64}())

    @test repr{Groebner.PowerVector{UInt8}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt8}())
    @test repr{Groebner.PowerVector{UInt16}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt16}())
    @test repr{Groebner.PowerVector{UInt32}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt32}())
    @test repr{Groebner.PowerVector{UInt64}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.NotPacked{UInt64}())
    @test repr{Groebner.PackedPair1{T, UInt8}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt8}())
    @test repr{Groebner.PackedPair1{T, UInt16}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt16}())
    @test repr{Groebner.PackedPair2{T, UInt32}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt32}())
    @test repr{Groebner.PackedPair3{T, UInt64}}() == Groebner.guess_effective_representation(f, unsaferepr(), Groebner.Lex(), Groebner.Packed{UInt64}())
end

@testset "f4 monom representations" begin
    using AbstractAlgebra
    for EV in all_hint_types
        for B in all_hint_params
            for domain in (GF(2^31-1), QQ)
                for s in [
                    Groebner.cyclicn(2, ground=domain),
                    Groebner.change_ordering(Groebner.noonn(4, ground=domain), :degrevlex),
                    Groebner.change_ordering(Groebner.katsuran(5, ground=domain), :degrevlex),
                    # Groebner.change_ordering(Groebner.s9_1(ground=domain), :degrevlex),
                    Groebner.change_ordering(Groebner.kinema(ground=domain), :degrevlex)
                ]
                    gb = Groebner.groebner(s, monoms=EV{B}())
                    @test Groebner.isgroebner(gb)
                end
            end
        end
    end
end
