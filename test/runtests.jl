
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

    @includetests ["core_f4_reduce", "core_f4_stress",
                    "core_f4", "rational_f4", "groebner_certify"]

    @includetests ["core_isgroebner", "core_isgroebner_stress",
                    "isgroebner_certify"]

    @includetests ["core_normalform", "core_normalform_stress"]

    @includetests ["hard_problems_f4", "large_problems_f4"]

    @includetests ["probabilistic_linalg"]

    if try_import(:DynamicPolynomials)
        @includetests ["io_dynamic", "core_dynamic"]
    end

    @includetests ["onthefly_order_change"]

    @includetests ["core_fglm", "core_kbase"]

    @includetests ["benchmark_systems"]

end
