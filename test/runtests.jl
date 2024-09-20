using Test, TestSetExtensions
using InteractiveUtils, Random

using AbstractAlgebra
using Groebner

# TODO: test examples in README.md (https://github.com/thchr/TestReadme.jl)
# TODO: test examples in the documentation

# Check invariants during testing.
Groebner.invariants_enabled() = true

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
    @time @includetests ["arithmetic"]

    # Different implementations of a monomial 
    @time @includetests ["monoms/exponentvector", "monoms/packedtuples"]
    @time @includetests ["monoms/monom_arithmetic", "monoms/monom_orders"]

    @time @includetests ["groebner"]

    @time @includetests ["learn_and_apply"]

    @time @includetests ["isgroebner"]

    @time @includetests ["normalform"]

    # Test for different frontends: 
    # - AbstractAlgebra.jl
    # - Nemo.jl
    # - DynamicPolynomials.jl
    @time @includetests ["input_output/AbstractAlgebra"]
    if try_import(:DynamicPolynomials)
        @time @includetests ["input_output/GroebnerDynamicPolynomialsExt"]
    end
    if try_import(:Nemo)
        @time @includetests ["input_output/Nemo"]
    end

    @time @includetests ["output_inferred"]

    # test for regressions
    @time @includetests ["regressions"]
end
