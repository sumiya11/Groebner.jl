using Test
using TestSetExtensions

using AbstractAlgebra
using Random
using Groebner

# Check invariants during testing.
# NOTE: it's good to turn this on!
Groebner.invariants_enabled() = true
Groebner.logger_update(loglevel=0)

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

# @test isempty(Test.detect_unbound_args(Groebner))
@test isempty(Test.detect_ambiguities(Groebner))

⊂(xs, ys) = all(in(ys), xs)
≂(xs, ys) = ⊂(xs, ys) && ⊂(ys, xs)

# Groebner.versioninfo()

@time @testset "All tests" verbose = true begin
    # Different implementations of a monomial 
    @includetests ["monoms/exponentvector", "monoms/packedtuples", "monoms/sparsevector"]
    # High-level monomial arithmetic and term orders
    @includetests ["monoms/monom_arithmetic", "monoms/monom_orders"]

    # Basic tests for addition in Zp
    @includetests ["arithmetic/Zp"]

    # Consistency of input-output
    @includetests ["input-output/AbstractAlgebra"]
    # Crt and rational reconstructions
    @includetests [
        "reconstruction/crt_reconstruction",
        "reconstruction/rational_reconstruction"
    ]

    @includetests [
        "groebner/groebner",
        "groebner/groebner_stress",
        "groebner/groebner_large",
        "groebner/many_variables",
        "groebner/large_exponents",
        "groebner/homogenization",
        "groebner/multi_threading"
    ]

    @includetests ["learn_and_apply/learn_and_apply", "learn_and_apply/apply_in_batches"]

    @includetests ["isgroebner/isgroebner"]

    @includetests ["normalform/normalform", "normalform/normalform_stress"]
    @includetests ["fglm/kbase"]

    # Test for different frontends: 
    # - AbstractAlgebra.jl  (AbstractAlgebra.Generic.MPoly{T})
    # - Nemo.jl  (Nemo.fmpq_mpoly, Nemo.gfp_mpoly)
    # - DynamicPolynomials.jl  (DynamicPolynomials.Polynomial{true, T})
    if try_import(:DynamicPolynomials)
        @includetests ["input-output/DynamicPolynomials"]
    end
    if try_import(:Nemo)
        @includetests ["input-output/Nemo"]
    end

    @includetests ["output_inferred"]

    @includetests ["utils/logging", "utils/timings"]

    # test for regressions
    @includetests ["regressions"]
end
