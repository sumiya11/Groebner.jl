
using Test
using TestSetExtensions

include("../src/GroebnerBases.jl")

using .GroebnerBases
using .GroebnerBases.AbstractAlgebra

@info "Testing started"

@testset "All tests" begin

    @includetests ARGS

end
