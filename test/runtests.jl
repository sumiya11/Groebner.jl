using Test
using InteractiveUtils, Random

using AbstractAlgebra
using Groebner

using AbstractAlgebra, Groebner

R, (a, b, beta, c, d, h, k, lm, q, u, v₁, v₂, v₃, v₄, v_₁, v_₂, v_₃, v_₄, v__₁, v__₂, v__₃, v__₄, w_₁, w_₂, w_₃, w_₄, w__₁, w__₂, w__₃, w__₄, x_₁, x_₂, x_₃, x_₄, x__₁, x__₂, x__₃, x__₄, y₁, y₂, y₃, y₄, y_₁, y_₂, y_₃, y_₄, y__₁, y__₂, y__₃, y__₄, z_₁, z_₂, z_₃, z_₄, z__₁, z__₂, z__₃, z__₄)=polynomial_ring(QQ,["a"," b"," beta"," c"," d"," h"," k"," lm"," q"," u"," v₁"," v₂"," v₃"," v₄"," v_₁"," v_₂"," v_₃"," v_₄"," v__₁"," v__₂"," v__₃"," v__₄"," w_₁"," w_₂"," w_₃"," w_₄"," w__₁"," w__₂"," w__₃"," w__₄"," x_₁"," x_₂"," x_₃"," x_₄"," x__₁"," x__₂"," x__₃"," x__₄"," y₁"," y₂"," y₃"," y₄"," y_₁"," y_₂"," y_₃"," y_₄"," y__₁"," y__₂"," y__₃"," y__₄"," z_₁"," z_₂"," z_₃"," z_₄"," z__₁"," z__₂"," z__₃"," z__₄"])

sys=[40-53*x__₁-1060*x_₁, 20-27*x__₂-540*x_₂, 100-137*x__₃-2740*x_₃, 200*y₂-y__₁-20*y_₁-200*y₁, 200*y₃-y__₂-20*y_₂-200*y₂, 200*y₄-y__₃-20*y_₃-200*y₃, 200*v₂-v__₁-20*v_₁-200*v₁, 200*v₃-v__₂-20*v_₂-200*v₂, 200*v₄-v__₃-20*v_₃-200*v₃, -25-2*w__₁-40*w_₁, -100-9*w__₂-180*w_₂, -10-w__₃-20*w_₃, -40-3*z__₁-60*z_₁, -200-17*z__₂-340*z_₂, -100-9*z__₃-180*z_₃, beta*v₁+d-6*lm+6*x_₁, 6*beta*v₁*x_₁+beta*v_₁+6*d*x_₁+6*x__₁, beta*v₂+d-6*lm+6*x_₂, 6*beta*v₂*x_₂+beta*v_₂+6*d*x_₂+6*x__₂, beta*v₃+d-6*lm+6*x_₃, 6*beta*v₃*x_₃+beta*v_₃+6*d*x_₃+6*x__₃, beta*v₄+d-6*lm+6*x_₄, 6*beta*v₄*x_₄+beta*v_₄+6*d*x_₄+6*x__₄, 6*a*y₁-beta*v₁+6*y_₁, -6*beta*v₁*x_₁+6*a*y_₁-beta*v_₁+6*y__₁, 6*a*y₂-beta*v₂+6*y_₂, -6*beta*v₂*x_₂+6*a*y_₂-beta*v_₂+6*y__₂, 6*a*y₃-beta*v₃+6*y_₃, -6*beta*v₃*x_₃+6*a*y_₃-beta*v_₃+6*y__₃, 6*a*y₄-beta*v₄+6*y_₄, -6*beta*v₄*x_₄+6*a*y_₄-beta*v_₄+6*y__₄, -k*y₁+u*v₁+v_₁, -k*y_₁+u*v_₁+v__₁, -k*y₂+u*v₂+v_₂, -k*y_₂+u*v_₂+v__₂, -k*y₃+u*v₃+v_₃, -k*y_₃+u*v_₃+v__₃, -k*y₄+u*v₄+v_₄, -k*y_₄+u*v_₄+v__₄, 6*c*q*y₁-c*y₁+6*b+9*w_₁, 18*c*q*w_₁*y₁+12*c*q*y_₁-3*c*w_₁*y₁-12*c*x_₁*y₁+18*b*w_₁-2*c*y_₁+18*w__₁, 6*c*q*y₂-c*y₂+6*b+9*w_₂, 18*c*q*w_₂*y₂+12*c*q*y_₂-3*c*w_₂*y₂-12*c*x_₂*y₂+18*b*w_₂-2*c*y_₂+18*w__₂, 50*c*q*y₃-9*c*y₃+50*b+90*w_₃, 90*c*q*w_₃*y₃+50*c*q*y_₃-15*c*w_₃*y₃-50*c*x_₃*y₃+90*b*w_₃-9*c*y_₃+90*w__₃, 11*c*q*y₄-2*c*y₄+11*b+22*w_₄, 66*c*q*w_₄*y₄+33*c*q*y_₄-11*c*w_₄*y₄-33*c*x_₄*y₄+66*b*w_₄-6*c*y_₄+66*w__₄, -10*c*q*y₁+12*h+15*z_₁, -3*c*q*w_₁*y₁-2*c*q*y_₁+3*h*z_₁+3*z__₁, -8*c*q*y₂+9*h+12*z_₂, -3*c*q*w_₂*y₂-2*c*q*y_₂+3*h*z_₂+3*z__₂, -5*c*q*y₃+6*h+9*z_₃, -9*c*q*w_₃*y₃-5*c*q*y_₃+9*h*z_₃+9*z__₃, -3*c*q*y₄+4*h+6*z_₄, -2*c*q*w_₄*y₄-c*q*y_₄+2*h*z_₄+2*z__₄, 4-5*v₁-5*y₁, 4-5*v₂-5*y₂, 3-4*v₃-4]

@info "" sys
@time groebner(sys)


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
