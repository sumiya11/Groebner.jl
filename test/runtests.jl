using Test
using TestSetExtensions

include("../src/Groebner.jl")

using .Groebner
using .Groebner.AbstractAlgebra

# if some particular time consuming tests should be executed
long_tests() = false

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

@time @testset "All tests" verbose=true begin
    @warn "Warnings during testing are fine."

    # test reference implementation
    @includetests ["reference"]
    # test different implementations of monomial 
    @includetests ["monoms/powervector", "monoms/packedpairs"]
    # test high-level monomial arithmetic and term orders
    @includetests ["monoms/term_arithmetic", "monoms/term_orders"]
    # test consistency of input-output
    @includetests ["input-output/io"]
    # test univriate input
    @includetests ["input-output/univariate"]
    # test crt and rational reconstructions
    @includetests ["crt_reconstruction", "mod_reconstruction"]
    # test `groebner`
    @includetests ["f4_reduce", 
                 "f4_stress", "adaptive_coefficients",
                 "f4", "rational_f4", "groebner_certify"]
    # test additional options in `groebner`
    @includetests ["probabilistic_linalg", "onthefly_order_change",
                "monom_representations"]
    # test some problems with very large output
    @includetests ["large_problems_f4"]
    # test `isgroebner`
    @includetests ["isgroebner", "isgroebner_stress", "isgroebner_certify"]
    # test `normalform` 
    @includetests ["normalform", "normalform_stress", "array_normalform"]
    # test `fglm` and `kbase`
    @includetests ["fglm", "kbase"]
    # test some special cases
    @includetests ["handling_zeros", "handling_checks", 
                "many_variables", "large_exponents"]
    if try_import(:DynamicPolynomials)
        @includetests ["input-output/io_dynamic", "dynamic"]
    end
    if try_import(:Nemo)
        @includetests ["input-output/io_nemo"]
    end

    # test for regressions
    @includetests ["regressions"]

    # test for some systems used in benchmarks
    @includetests ["small_benchmarks"]
    if long_tests()
        @includetests ["large_benchmarks"]
    end
end
