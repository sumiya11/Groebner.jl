
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

    @includetests ["commons", "io", "io_consistency",
                    "univariate_aa"]

    @includetests ["crt_reconstruction"] # "mod_reconstruction"]#,

    @includetests ["f4_reduce", "f4_stress",
                    "f4", "rational_f4", "groebner_certify"]

    @includetests ["isgroebner", "isgroebner_stress",
                    "isgroebner_certify"]

    @includetests ["normalform", "normalform_stress"]

    @includetests ["hard_problems_f4", "large_problems_f4"]

    @includetests ["probabilistic_linalg"]

    if try_import(:DynamicPolynomials)
        @includetests ["io_dynamic", "dynamic"]
    end

    @includetests ["onthefly_order_change"]

    @includetests ["fglm", "kbase"]

    @includetests ["benchmark_systems"]

end
