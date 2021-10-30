
using .GroebnerBases: f4, rootn, katsura6

@testset "F4 over Finite Fields" begin

    # Root n
    ground = GF(2^31-1)
    fs = rootn(3, ground=ground)
    x1, x2, x3 = gens(parent(first(fs)))

    G = f4(fs, reduced=true)

    @test G == [x3^3 - 1, x2^2 + x2*x3 + x3^2, x1 + x2 + x3]

    # katsura 6
    # a
    
end
