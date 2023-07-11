using Test
using TestSetExtensions

using AbstractAlgebra
using Groebner

# Check invariants during testing
Groebner.invariants_enabled() = true

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

@test isempty(Test.detect_unbound_args(Groebner))
@test isempty(Test.detect_ambiguities(Groebner))

⊂(xs, ys) = all(in(ys), xs)
≂(xs, ys) = ⊂(xs, ys) && ⊂(ys, xs)

@time @testset "All tests" verbose = true begin
    # Different implementations of a monomial 
    @includetests ["monoms/exponentvector", "monoms/packedtuples"]
    # High-level monomial arithmetic and term orders
    @includetests ["monoms/monom_arithmetic", "monoms/monom_orders"]

    # Consistency of input-output
    @includetests ["input-output/AbstractAlgebra"]
    # Crt and rational reconstructions
    @includetests ["groebner/crt_reconstruction", "groebner/mod_reconstruction"]

    @includetests [
        "groebner/groebner",
        "groebner/groebner_stress",
        "groebner/groebner_large",
        "groebner/many_variables",
        "groebner/large_exponents"
    ]

    @includetests ["isgroebner/isgroebner"]

    @includetests ["normalform/normalform", "normalform/normalform_stress"]
    @includetests ["fglm/kbase"]

    # Groebner.jl is tested for different frontends: 
    # - AbstractAlgebra.jl  (AbstractAlgebra.Generic.MPoly{T})
    # - Nemo.jl  (Nemo.fmpq_mpoly, Nemo.gfp_mpoly, ...)
    # - DynamicPolynomials.jl (DynamicPolynomials.Polynomial{true, T})
    if try_import(:DynamicPolynomials)
        @includetests ["input-output/DynamicPolynomials"]
    end
    if try_import(:Nemo)
        @includetests ["input-output/Nemo"]
    end

    @includetests ["output_inferred"]

    # test for regressions
    @includetests ["regressions"]
end
