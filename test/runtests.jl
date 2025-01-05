using Test
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
    @time include("arithmetic.jl")

    # Different implementations of a monomial 
    @time include("monoms/exponentvector.jl")
    @time include("monoms/packedtuples.jl")
    @time include("monoms/monom_arithmetic.jl")
    @time include("monoms/monom_orders.jl")

    @time include("groebner.jl")

    @time include("learn_and_apply.jl")

    @time include("isgroebner.jl")

    @time include("normalform.jl")

    @time include("auxiliary.jl")

    # Test for different frontends: 
    # - AbstractAlgebra.jl
    # - Nemo.jl
    # - DynamicPolynomials.jl
    @time include("input_output/AbstractAlgebra.jl")
    if try_import(:DynamicPolynomials)
        @info "Testing frontend: DynamicPolynomials.jl"
        @time include("input_output/GroebnerDynamicPolynomialsExt.jl")
    end
    if try_import(:Nemo)
        @info "Testing frontend: Nemo.jl"
        @time include("input_output/Nemo.jl")
    end

    @time include("output_inferred.jl")

    # test for regressions
    @time include("regressions.jl")
end
