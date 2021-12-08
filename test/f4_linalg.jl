

using .GroebnerBases: constructrows_sparse, rref_sparse!,
                    constructpolys_sparse, linear_algebra_sparse_det,
                    linear_algebra_dense_det, constructpolys_dense,
                    constructrows_dense

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

#=
@testset "F4 construct matrix" begin

    @test constructmatrix([R(4)], [R(1)]) ==  reshape([4], 1, 1)
    @test constructmatrix([x], [x]) == reshape([1], 1, 1)
    @test constructmatrix([5x + y], [x, y]) == [5 1]

    fs = [x^2 - 3y^2, x + 3y^2 + 1, (4//3)z^3]
    Ts = sort(monomials(fs), rev=true)
    @test constructmatrix(fs, Ts) == [
        1 0 -3 0 0;
        0 1  3 0 1;
        0 0  0 QQ(4, 3) 0;
    ]

    fs = [x^3 + x^2, x^3 - x^2, x^3, x^2]
    Ts = sort(monomials(fs), rev=true)
    @test constructmatrix(fs, Ts) == [
        1  1;
        1 -1;
        1  0;
        0  1;
    ]

end

@testset "F4 construct polys" begin

    fs = [5x]
    Ts = sort(monomials(fs), rev=true)
    A = reshape(QQ.([5]), (1, 1))
    @test constructpolys((1, A), Ts) == fs

    fs = [x^2 - y + 4, x + 2, 2y, R(8)]
    Ts = sort(monomials(fs), rev=true)
    A = QQ.([
        1 0 -1 4;
        0 1  0 2;
        0 0  2 0;
        0 0 0  8;
    ])
    @test constructpolys((4, A), Ts) == fs

    fs = [x^3 + x^2, -x^2]
    Ts = sort(monomials(fs), rev=true)
    A = QQ.([
        1  1;
        0 -1;
        0  0;
        0  0;
    ])
    @test constructpolys((2, A), Ts) ==  fs

    # check that polys -> matrix -> polys works
    fs = [x^5 + x - y + 1, z^2 - 2 + x^4,
        x^4 - x^3 + x^2 + 1, x^2 + y^2 + z^2]
    Ts = sort(monomials(fs), rev=true)
    @test constructpolys((4, constructmatrix(fs, Ts)), Ts) == fs

    # check that matrix -> polys -> matrix works
    A = QQ.([
        1 2 3 4 0 0;
        4 2 0 1 1 -1;
        0 0 0 0 0 1;
        0 1 1 0 1 0;
        1 1 1 1 1 1;
        0 0 0 2 0 0;
    ])
    Ts = sort([x^3, x, y^2, z^2, z, R(1)], rev=true)
    @test constructmatrix(constructpolys((6, A), Ts), Ts) == A

end

@testset "F4 linear algebra reduction" begin
    @test linear_algebra_aa([x]) ⊂ [x]
    @test linear_algebra_aa([x, x^2]) ⊂ [x^2, x]
    @test linear_algebra_aa([x + 1, x, R(1)]) ⊂ [x, 1]

    @test linear_algebra_aa([x + 5, x - 5y]) ⊂ [x + 5, y + 1]
    @test linear_algebra_aa([x + y + z, y + z, z]) ⊂ [x, y, z]

    fs = [x^2 + y^2 + z^2 + 1, x^2 + 2y^2, R(1)]
    @test linear_algebra_aa(fs) ⊂ [y^2 - z^2, x^2 + 2z^2, R(1)]

end
=#

#=
@testset "F4 construct sparse rows" begin

    @test constructrows_sparse([R(4)], [R(1)]) == [sparse([QQ(4)])]
    @test constructrows_sparse([x], [x]) == [sparse([QQ(1)])]
    @test constructrows_sparse([5x + y], [x, y]) == [sparse([QQ(5), QQ(1)])]

    fs = [x^2 - 3y^2, x + 3y^2 + 1, (4//3)z^3]
    Ts = sort(monomials(fs), rev=true)
    @test constructrows_sparse(fs, Ts) == [
        sparse(QQ.([0, 0,  0, 4//3, 0])),
        sparse(QQ.([0, 1,  3, 0, 1])),
        sparse(QQ.([1, 0, -3, 0, 0])),
    ]

    fs = [x^3 + x^2, x^3 - x^2, x^3, x^2]
    Ts = sort(monomials(fs), rev=true)
    @test constructrows_sparse(fs, Ts) == [
        sparse([0,  1]),
        sparse([1,  0]),
        sparse([1,  1]),
        sparse([1, -1]),
    ]

end
=#

#=
@testset "F4 construct polys from sparse" begin

    fs = [5x]
    Ts = sort(monomials(fs), rev=true)
    A = reshape(QQ.([5]), (1, 1))
    @test constructpolys((1, A), Ts) == fs

    fs = [x^2 - y + 4, x + 2, 2y, R(8)]
    Ts = sort(monomials(fs), rev=true)
    A = QQ.([
        1 0 -1 4;
        0 1  0 2;
        0 0  2 0;
        0 0 0  8;
    ])
    @test constructpolys((4, A), Ts) == fs

    fs = [x^3 + x^2, -x^2]
    Ts = sort(monomials(fs), rev=true)
    A = QQ.([
        1  1;
        0 -1;
        0  0;
        0  0;
    ])
    @test constructpolys((2, A), Ts) ==  fs

    # check that polys -> matrix -> polys works
    fs = [x^5 + x - y + 1, z^2 - 2 + x^4,
        x^4 - x^3 + x^2 + 1, x^2 + y^2 + z^2]
    Ts = sort(monomials(fs), rev=true)
    @test constructpolys((4, constructmatrix(fs, Ts)), Ts) == fs

    # check that matrix -> polys -> matrix works
    A = QQ.([
        1 2 3 4 0 0;
        4 2 0 1 1 -1;
        0 0 0 0 0 1;
        0 1 1 0 1 0;
        1 1 1 1 1 1;
        0 0 0 2 0 0;
    ])
    Ts = sort([x^3, x, y^2, z^2, z, R(1)], rev=true)
    @test constructmatrix(constructpolys((6, A), Ts), Ts) == A

end
=#

#=
@testset "F4 linear algebra sparse det" begin
    @test GroebnerBases.linear_algebra_sparse_det([x]) ⊂ [x]
    @test GroebnerBases.linear_algebra_sparse_det([x, x^2]) ⊂ [x^2, x]
    @test GroebnerBases.linear_algebra_sparse_det([x + 1, x, R(1)]) ⊂ [x, 1]

    @test GroebnerBases.linear_algebra_sparse_det([x + 5, x - 5y]) ⊂ [x + 5, y + 1]
    @test GroebnerBases.linear_algebra_sparse_det([x + y + z, y + z, z]) ⊂ [x, y, z]

    fs = [x^2 + y^2 + z^2 + 1, x^2 + 2y^2, R(1)]
    # WE DO NOT PERFORM BACKWARD SUBSTITUTION
    @test GroebnerBases.linear_algebra_sparse_det(fs) ⊂ [y^2 - z^2, x^2 + 2y^2, R(1)]

end
=#

#=
@testset "F4 linear algebra dense det" begin
    @test GroebnerBases.linear_algebra_dense_det([x]) ⊂ [x]
    @test GroebnerBases.linear_algebra_dense_det([x, x^2]) ⊂ [x^2, x]
    @test GroebnerBases.linear_algebra_dense_det([x + 1, x, R(1)]) ⊂ [x, 1]

    @test GroebnerBases.linear_algebra_dense_det([x + 5, x - 5y]) ⊂ [x + 5, y + 1]
    @test GroebnerBases.linear_algebra_dense_det([x + y + z, y + z, z]) ⊂ [x, y, z]

    fs = [x^2 + y^2 + z^2 + 1, x^2 + 2y^2, R(1)]
    # WE DO NOT PERFORM BACKWARD SUBSTITUTION
    @test GroebnerBases.linear_algebra_dense_det(fs) ⊂ [y^2 - z^2, x^2 + 2y^2, R(1)]
end
=#

#=
@testset "F4 linear algebra reduction" begin
    @test linear_algebra_aa([x]) ⊂ [x]
    @test linear_algebra_aa([x, x^2]) ⊂ [x^2, x]
    @test linear_algebra_aa([x + 1, x, R(1)]) ⊂ [x, 1]

    @test linear_algebra_aa([x + 5, x - 5y]) ⊂ [x + 5, y + 1]
    @test linear_algebra_aa([x + y + z, y + z, z]) ⊂ [x, y, z]

    fs = [x^2 + y^2 + z^2 + 1, x^2 + 2y^2, R(1)]
    @test linear_algebra_aa(fs) ⊂ [y^2 - z^2, x^2 + 2z^2, R(1)]

end
=#
