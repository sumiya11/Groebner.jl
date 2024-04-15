# This script is not intended to be run directly.
# Instead, see runtests.jl.

function dump_results(file, key)
    open(file, "w") do out
        println(out, key)
        for problem in suite
            problem_name = problem.problem_name
            result = problem.result
            println(out, "$problem_name: $(join(map(string, result), ","))")
        end
    end
end

# TTFX
t1 = @timed using Groebner
using AbstractAlgebra

R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering=:degrevlex)
sys = [x^2 + y, x * y^2]
t2 = @timed groebner(sys)

R, (x, y) = polynomial_ring(GF(2^30 + 3), ["x", "y"])
sys = [x^2 + y, x * y^2]
t3 = @timed groebner(sys)

suite = []

push!(suite, (problem_name="using Groebner", result=[t1.time]))
push!(suite, (problem_name="groebner, ttfx, qq drl", result=[t2.time]))
push!(suite, (problem_name="groebner, ttfx, zp lex", result=[t3.time]))

# Allocations
k = Groebner.katsuran(3, internal_ordering=:degrevlex)
groebner(k)
a1 = @allocated groebner(k)

k = Groebner.katsuran(8, k=GF(2^30 + 3), internal_ordering=:degrevlex)
groebner(k)
a2 = @allocated groebner(k)

# push!(suite, (problem_name="groebner, allocs qq", result=[a1]))
# push!(suite, (problem_name="groebner, allocs zp", result=[a2]))

# Runtime
import Primes
import Nemo

function nemo_make_prime_finite_field(p)
    if p < typemax(UInt)
        Nemo.fpField(convert(UInt, p), false)
    else
        Nemo.FpField(Nemo.ZZRingElem(p), false)
    end
end

function compute_gb(system, trials=7; kws...)
    times = []
    for _ in 1:trials
        GC.gc()
        time = @elapsed groebner(system; kws...)
        push!(times, time)
    end
    times
end

# Compute Groebner bases over integers modulo a large prime
problem = (
    problem_name="groebner, AA, GF(2^31-1), katsura 5",
    result=compute_gb(Groebner.katsuran(5, internal_ordering=:degrevlex, k=GF(2^31 - 1)))
)
push!(suite, problem)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 5",
        result=compute_gb(
            Groebner.katsuran(5, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 6",
        result=compute_gb(
            Groebner.katsuran(6, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^31-1), katsura 8",
#         result=compute_gb(
#             Groebner.katsuran(8, internal_ordering=:degrevlex, k=GF(2^31 - 1))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^31-1), katsura 10",
#         result=compute_gb(
#             Groebner.katsuran(10, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
#             5
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^27+29), katsura 10",
#         result=compute_gb(
#             Groebner.katsuran(10, internal_ordering=:degrevlex, k=GF(2^27 + 29)),
#             5
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^30+3), katsura 11",
#         result=compute_gb(
#             Groebner.katsuran(11, internal_ordering=:degrevlex, k=GF(2^30 + 3)),
#             3
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^27+29), cyclic 8",
#         result=compute_gb(
#             Groebner.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^27 + 29)),
#             5
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^31-1), cyclic 8",
#         result=compute_gb(
#             Groebner.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
#             5
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, Nemo, GF(2^31-1), cyclic 8",
#         result=compute_gb(
#             Groebner.cyclicn(
#                 8,
#                 internal_ordering=:degrevlex,
#                 k=nemo_make_prime_finite_field(2^31 - 1)
#             )
#         )
#     )
# )

# function n_variable_set(n; internal_ordering=:degrevlex, k=GF(2^31 - 1))
#     R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
#     f = [sum(prod(x[i:(n - kk)], init=1) for i in 1:(kk + 1)) for kk in 0:(n - 1)]
#     f
# end

# push!(
#     suite,
#     (
#         problem_name="groebner, AA, GF(2^31-1), 100 vars",
#         result=compute_gb(
#             n_variable_set(100, internal_ordering=:degrevlex, k=GF(2^31 - 1))
#         )
#     )
# )

# push!(
#     suite,
#     (
#         problem_name="groebner, threaded, AA, GF(2^31-1), cyclic 8",
#         result=compute_gb(
#             Groebner.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
#             5,
#             threaded=:yes
#         )
#     )
# )

# function learn_and_apply(system)
#     times = []
#     trials = 7
#     graph, gb = groebner_learn(system)
#     for _ in 1:trials
#         GC.gc()
#         time2 = @elapsed groebner_apply!(graph, system)
#         push!(times, time2)
#     end
#     times
# end

# # Learn and apply
# push!(
#     suite,
#     (
#         problem_name="groebner_apply!, AA, GF(2^31-1), cyclic 7",
#         result=learn_and_apply(
#             Groebner.cyclicn(7, internal_ordering=:degrevlex, k=GF(2^31 - 1))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner_apply!, Nemo, GF(2^31-1), cyclic 7",
#         result=learn_and_apply(
#             Groebner.cyclicn(
#                 7,
#                 internal_ordering=:degrevlex,
#                 k=nemo_make_prime_finite_field(2^31 - 1)
#             )
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner_apply!, AA, GF(2^31-1), katsura 10",
#         result=learn_and_apply(
#             Groebner.katsuran(10, internal_ordering=:degrevlex, k=GF(2^31 - 1))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner_apply!, AA, GF(2^27+29), katsura 10",
#         result=learn_and_apply(
#             Groebner.katsuran(10, internal_ordering=:degrevlex, k=GF(2^27 + 29))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner_apply!, Nemo, GF(2^31-1), katsura 10",
#         result=learn_and_apply(
#             Groebner.katsuran(
#                 10,
#                 internal_ordering=:degrevlex,
#                 k=nemo_make_prime_finite_field(2^31 - 1)
#             )
#         )
#     )
# )

# # Compute Groebner bases over the rationals
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, katsura 8",
#         result=compute_gb(Groebner.katsuran(8, internal_ordering=:degrevlex, k=QQ))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, Nemo, QQ, katsura 8",
#         result=compute_gb(Groebner.katsuran(8, internal_ordering=:degrevlex, k=Nemo.QQ))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, eco 10",
#         result=compute_gb(Groebner.eco10(internal_ordering=:degrevlex, k=QQ))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, cyclic 7",
#         result=compute_gb(Groebner.cyclicn(7, internal_ordering=:degrevlex, k=QQ))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, noon 8",
#         result=compute_gb(Groebner.noonn(8, internal_ordering=:degrevlex, k=QQ), 3)
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, hexapod",
#         result=compute_gb(Groebner.hexapod(internal_ordering=:degrevlex, k=QQ), 3)
#     )
# )

# function multimodular_gb_problem(nbits; np=AbstractAlgebra)
#     R, (x1, x2, x3, x4) =
#         polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)
#     nbits_per_prime = 31
#     nprimes = max(div(nbits, nbits_per_prime), 1)
#     N = prod(map(BigInt, Primes.nextprimes(2^31 - 100, nprimes)))
#     @info "Constructing a multi-modular problem" np nbits nbits_per_prime nprimes
#     system = [
#         x1 + x2 + x3 + x4,
#         x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
#         x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
#         x1 * x2 * x3 * x4 + N
#     ]
#     system
# end

# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, 1000-bit output coeffs",
#         result=compute_gb(multimodular_gb_problem(1000, np=AbstractAlgebra))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, AA, QQ, 10000-bit output coeffs",
#         result=compute_gb(multimodular_gb_problem(10000, np=AbstractAlgebra))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, Nemo, QQ, 1000-bit output coeffs",
#         result=compute_gb(multimodular_gb_problem(1000, np=Nemo))
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="groebner, Nemo, QQ, 10000-bit output coeffs",
#         result=compute_gb(multimodular_gb_problem(10000, np=Nemo))
#     )
# )

# function compute_normalforms(system)
#     R = parent(system[1])
#     gb = Groebner.groebner(system)
#     times = []
#     trials = 7
#     for _ in 1:trials
#         GC.gc()
#         time = @elapsed begin
#             n1 = normalform(gb, system)
#             n2 = normalform(gb, gb)
#         end
#         push!(times, time)
#     end
#     times
# end

# # Normal forms
# push!(
#     suite,
#     (
#         problem_name="normalform, AA, GF(2^31-1), cyclic 7",
#         result=compute_normalforms(
#             Groebner.cyclicn(7, internal_ordering=:degrevlex, k=GF(2^31 - 1))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="normalform, AA, GF(103), cyclic 8",
#         result=compute_normalforms(
#             Groebner.cyclicn(8, internal_ordering=:degrevlex, k=GF(103))
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="normalform, Nemo, GF(103), cyclic 8",
#         result=compute_normalforms(
#             Groebner.cyclicn(
#                 8,
#                 internal_ordering=:degrevlex,
#                 k=nemo_make_prime_finite_field(103)
#             )
#         )
#     )
# )
# push!(
#     suite,
#     (
#         problem_name="normalform, AA, QQ, katsura 9",
#         result=compute_normalforms(
#             Groebner.katsuran(9, internal_ordering=:degrevlex, k=QQ)
#         )
#     )
# )
