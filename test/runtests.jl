
using Test
using TestSetExtensions

include("../src/GroebnerBases.jl")

using .GroebnerBases
using .GroebnerBases.AbstractAlgebra

# @info "Testing started"

⊂(xs, ys) = all(in(ys), xs)
# set equality
≂(xs, ys) = ⊂(xs, ys) && ⊂(ys, xs)

@testset "All tests" begin

    @includetests
end
