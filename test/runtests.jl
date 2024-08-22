using Test, TestSetExtensions
using InteractiveUtils, Random

using AbstractAlgebra
using Groebner

# TODO: test examples in README.md (https://github.com/thchr/TestReadme.jl)
# TODO: test examples in the documentation

# Check invariants during testing.
Groebner.invariants_enabled() = true
Groebner.logging_enabled() = true

function is_github_ci()
    return parse(Bool, get(ENV, "GITHUB_ACTIONS", "false"))
end

if is_github_ci()
    @info "Running in a Github CI job. Printing specs!"
    versioninfo(verbose=true)
else
    versioninfo()
end

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

@test isempty(Test.detect_ambiguities(Groebner))

@time @testset "All tests" verbose = true begin
    # Basic tests for addition in Zp
    @time @includetests ["arithmetic/Zp"]

    # Different implementations of a monomial 
    @time @includetests [
        "monoms/exponentvector",
        "monoms/packedtuples",
        "monoms/sparsevector"
    ]
    # High-level monomial arithmetic and term orders
    @time @includetests ["monoms/monom_arithmetic", "monoms/monom_orders"]

    # Consistency of input-output
    @time @includetests ["input_output/AbstractAlgebra"]

    @time @includetests [
        "groebner/groebner",
        "groebner/groebner_stress",
        "groebner/homogenization",
        "groebner/multi_threading"
    ]

    @time @includetests ["groebner/groebner_with_change_matrix"]

    @time @includetests [
        "learn_and_apply/learn_and_apply",
        "learn_and_apply/apply_in_batches"
    ]

    @time @includetests ["isgroebner/isgroebner"]

    @time @includetests ["normalform/normalform", "normalform/normalform_stress"]

    # Test for different frontends: 
    # - AbstractAlgebra.jl  (AbstractAlgebra.Generic.MPoly{T})
    # - Nemo.jl  (Nemo.fmpq_mpoly, Nemo.gfp_mpoly)
    # - DynamicPolynomials.jl  (DynamicPolynomials.Polynomial{true, T})
    if try_import(:DynamicPolynomials)
        # @time @includetests ["input_output/DynamicPolynomials"]
    end
    if try_import(:Nemo)
        @time @includetests ["input_output/Nemo"]
    end

    @time @includetests ["output_inferred"]

    @time @includetests ["utils/logging"]

    # test for regressions
    @time @includetests ["regressions"]
end
