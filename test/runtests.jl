using Test, InteractiveUtils, Groebner
using InteractiveUtils, Random

# Check invariants during testing.
@eval Groebner invariants_enabled() = true

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

# Isolate tests into separate modules so they don't interfere
function isolated_testset(name)
    @assert endswith(name, ".jl")
    @info "Running $name"
    mod = gensym(Symbol(chop(name, tail=2)))
    @time @eval module $mod
    include($name)
    end
end

@test isempty(Test.detect_ambiguities(Groebner))

@time @testset "All tests" verbose = true begin
    isolated_testset("arithmetic.jl")
    isolated_testset("crt.jl")
    isolated_testset("ratrec.jl")

    # Different implementations of a monomial
    isolated_testset("monoms/exponentvector.jl")
    isolated_testset("monoms/packedtuples.jl")
    isolated_testset("monoms/monom_arithmetic.jl")
    isolated_testset("monoms/monom_orders.jl")

    isolated_testset("groebner.jl")
    isolated_testset("learn_and_apply.jl")
    isolated_testset("isgroebner.jl")
    isolated_testset("normalform.jl")
    isolated_testset("auxiliary.jl")

    isolated_testset("output_inferred.jl")
    isolated_testset("regressions.jl")

    # Test for different frontends: 
    # - AbstractAlgebra.jl
    # - Nemo.jl
    # - DynamicPolynomials.jl
    # isolated_testset("input_output/AbstractAlgebra.jl")
    if try_import(:DynamicPolynomials)
        isolated_testset("input_output/GroebnerDynamicPolynomialsExt.jl")
    end
    if try_import(:Nemo)
        isolated_testset("input_output/Nemo.jl")
    end
end
