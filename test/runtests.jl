
using Test
using TestSetExtensions

include("../src/Groebner.jl")

using .Groebner
using .Groebner.AbstractAlgebra

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

⊂(xs, ys) = all(in(ys), xs)
# set equality
≂(xs, ys) = ⊂(xs, ys) && ⊂(ys, xs)


@testset "All tests" begin

    @includetests ["commons", "io", "io_consistency"]

    @includetests ["coeffs", "sorting",
                    "mod_reconstruction", "crt_reconstruction"]

    @includetests ["core_f4_reduce", "core_f4_stress",
                    "core_f4", "rational_f4"]

    @includetests ["core_isgroebner", "core_isgroebner_stress"]

    @includetests ["core_normalform", "core_normalform_stress"]

    @includetests ["hard_problems_f4", "large_problems_f4"]

    if try_import(:DynamicPolynomials)
        @includetests ["dynamic_io", "dynamic_core"]
    end
end
