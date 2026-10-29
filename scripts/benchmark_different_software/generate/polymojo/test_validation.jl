using Test
include("validate.jl")

@testset "Complete output comparison" begin
    expected = "x,y\n0\nx + y,\ny^2 - 1\n"
    @test compare_full_bases(expected, "x,y\n0\ny^2 - 1,\ny + x\n")
    @test !compare_full_bases(expected, "x,y\n0\nx - y,\ny^2 - 1\n")
    @test !compare_full_bases(expected, "x,y\n0\ny^2 - 1\n")
    @test !compare_full_bases(expected, "x,y\n7\nx + y,\ny^2 - 1\n")
    @test !compare_full_bases(expected, "y,x\n0\nx + y,\ny^2 - 1\n")
    @test !compare_full_bases(expected, "x,y\n0\n2*x + 2*y,\ny^2 - 1\n")
    @test compare_full_bases("x\n0\nx + 1/2\n", "x\n0\nx + 2/4\n")
end

@testset "Reference outputs up to nonzero scalar factors" begin
    expected = "x,y\n0\nx + y,\ny^2 - 1\n"
    @test compare_full_bases(expected, "x,y\n0\n-3*y^2 + 3,\n2*x + 2*y\n"; normalize_monic=true)
    @test !compare_full_bases(expected, "x,y\n0\n2*x - 2*y,\ny^2 - 1\n"; normalize_monic=true)
    @test !compare_full_bases(expected, "x,y\n0\ny^2 - 1\n"; normalize_monic=true)
    @test !compare_full_bases(expected, "x,y\n7\nx + y,\ny^2 - 1\n"; normalize_monic=true)
    @test !compare_full_bases(expected, "y,x\n0\nx + y,\ny^2 - 1\n"; normalize_monic=true)
    @test compare_full_bases("x\n7\nx + 1\n", "x\n7\n3*x + 3\n"; normalize_monic=true)
    @test compare_full_bases("x\n0\n0\n", "x\n0\n0\n"; normalize_monic=true)
    @test !compare_full_bases("x\n0\nx\n", "x\n0\nx,\n2*x\n"; normalize_monic=true)
end
