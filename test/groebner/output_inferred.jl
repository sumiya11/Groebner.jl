
@testset "Output type inferred" begin
    using AbstractAlgebra
    R, (x, y) = QQ["x","y"]
    @inferred groebner([x, y])
    @inferred groebner([x, y], ordering=DegRevLex())
    @inferred normalform([x, y], x + 1)
    @inferred normalform([x, y], [x + 1, y + 1, R(0)])
    @inferred normalform([x, y], [R(0), R(0), R(0)])
    @inferred isgroebner([x, y])
    @inferred kbase([x^2, y])

    R, (x, y) = GF(2^31-1)["x","y"]
    @inferred groebner([x, y])
    @inferred groebner([x, y], ordering=DegRevLex())
    @inferred normalform([x, y], x + 1)
    @inferred normalform([x, y], [x + 1, y + 1, R(0)])
    @inferred normalform([x, y], [R(0), R(0), R(0)])
    @inferred isgroebner([x, y])
    @inferred kbase([x^2, y])
end
