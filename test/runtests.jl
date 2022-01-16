
using Test
using TestSetExtensions

include("../src/Groebner.jl")

using .Groebner
using .Groebner.AbstractAlgebra

# @info "Testing started"

⊂(xs, ys) = all(in(ys), xs)
# set equality
≂(xs, ys) = ⊂(xs, ys) && ⊂(ys, xs)

@testset "All tests" begin

    @includetests
end
