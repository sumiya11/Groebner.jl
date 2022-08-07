
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

    @includetests ["internal_term_orders"]

    @includetests ["reference", "io", "io_consistency",
                    "univariate_aa"]

    @includetests ["crt_reconstruction"] # "mod_reconstruction"]#,

    @includetests ["f4_reduce", "f4_stress", "adaptive_coefficients",
                    "f4", "rational_f4", "groebner_certify"]

    @includetests ["isgroebner", "isgroebner_stress",
                    "isgroebner_certify"]

    @includetests ["normalform", "normalform_stress",
                    "array_normalform"]

    @includetests ["hard_problems_f4", "large_problems_f4"]

    @includetests ["probabilistic_linalg"]

    if try_import(:DynamicPolynomials)
        @includetests ["io_dynamic", "dynamic"]
    end

    if try_import(:Nemo)
        @includetests ["io_nemo"]
    end

    @includetests ["onthefly_order_change", "handling_zeros",
                    "handling_checks"]

    @includetests ["fglm", "kbase"]

    @includetests ["regressions"]

    @includetests ["benchmark_systems"]

end
