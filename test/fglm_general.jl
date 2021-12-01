
using .GroebnerBases: fglm, groebner, f4, change_ordering,
                        reducegb, noon3, eco5


# we assume here that f4 is correct
@testset "fglm over Finite Fields" begin

    ground = GF(2^31-1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"], ordering=:degrevlex)

    fs_deg = [
        x^2 - 1,
        y^2 - 1,
        z^2 - y^3 - x^3 - 28
    ]

    gb_deg = f4(fs_deg)
    gb_lex_fglm = fglm(gb_deg)

    fs_lex = change_ordering(fs_deg, :lex)
    gb_lex = f4(fs_lex)

    @test reducegb(gb_lex) == reducegb(gb_lex_fglm)


    fs_deg = [
        x^2 + 2y^2 - y - 2z,
        x^2 - 8y^2 + 10z - 1,
        x^2 - 7y*z
    ]

    gb_deg = f4(fs_deg)
    gb_lex_fglm = fglm(gb_deg)

    fs_lex = change_ordering(fs_deg, :lex)
    gb_lex = f4(fs_lex)

    @test reducegb(gb_lex) == reducegb(gb_lex_fglm)


    # here we check fglm for
    # 5 zero dimensional systems of the same cyclic structure
    for i in 1:4
        fs_lex = rootn(i, ground=GF(2^31 - 1))
        fs_deg = change_ordering(fs_lex, :degrevlex)

        gb_lex = f4(fs_lex)

        gb_deg = f4(fs_deg)
        gb_lex_fglm = fglm(gb_deg)

        @test reducegb(gb_lex) == reducegb(gb_lex_fglm)
    end

    # NOON !
    #=
    noon = noon3(ground=ground)
    noongb_lex = f4(noon)

    noon_deg = change_ordering(noon, :degrevlex)
    noongb_deg = f4(noon_deg)
    noongb_fglm = fglm(noongb_deg)

    @test reducegb(noongb_fglm) == reducegb(noongb_lex)
    =#


    # Mickey-mouse example to illustrate homotopy continuation.
    R, (x, y) = PolynomialRing(ground, ["x", "y"], ordering=:lex)
    fs_lex = [
        x^2 + 4*y^2 - 4,
        2*y^2 - x
    ]
    fs_deg = change_ordering(fs_lex, :degrevlex)

    gb_lex = f4(fs_lex)
    gb_deg = fglm(f4(fs_deg))
    @test reducegb(gb_lex) == reducegb(gb_deg)

    #=
    R, (x, y, a, b) = PolynomialRing(ground, ["x", "y", "a", "b"], ordering=:lex)
    fs_lex = [
        3*x^2-2*x-a,
        x^3-x^2-x*a+a-2*b-2,
        3*y^2-2*y-a,
        y^3-y^2-y*a-a+2
    ]
    fs_deg = change_ordering(fs_lex, :degrevlex)

    gb_lex = f4(fs_lex)
    gb_deg = fglm(f4(fs_deg))
    @test reducegb(gb_lex) == reducegb(gb_deg)

    =#

end
