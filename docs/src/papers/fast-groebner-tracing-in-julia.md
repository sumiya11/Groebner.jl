# Groebner.jl: Fast Gröbner Tracing in Julia

This page mirrors the code snippets from the paper _Groebner.jl: Fast Gröbner Tracing in Julia_ that appeared in the proceedings of the International Congress on Mathematical Software 2026.

## Setup

```julia
using Pkg; Pkg.add("Groebner"); Pkg.add("Nemo")
```

## Generic coefficient arithmetic

This snippet matches the paper example computing a basis over a rational function field:

```@example paper_generic_coefficients
using Groebner, Nemo

K, t = rational_function_field(QQ, "t")
R, (x, y, z) = K["x", "y", "z"]

F = [t * x^2 * y + y^3, x * y^2 + 3 // t * z]

G = groebner(F)
isgroebner(G)
```

## Learn and apply

The paper uses the Katsura-10 system to illustrate the tracing interface over two primes:

```@example paper_trace
using Groebner, Nemo

kat = Groebner.Examples.katsuran(10)

kat_zp1 = map(f -> change_base_ring(GF(2^30 + 3), f), kat)
kat_zp2 = map(f -> change_base_ring(GF(2^30 + 7), f), kat)

trace, _ = groebner_learn(kat_zp1)
success, gb = groebner_apply!(trace, kat_zp2)

success && isgroebner(gb)
```

The `success` flag is a consistency check against the learned trace. In practice it is useful, but it is not a formal proof that the result is a Groebner basis.

## Timing the apply phase

The paper benchmarks these calls with `BenchmarkTools`. In the documentation we keep this to a single timed run so the page remains practical to execute in CI.

```@example paper_trace
trace_timing, _ = groebner_learn(kat_zp1)

groebner_time = @elapsed groebner(kat_zp2)
apply_time = @elapsed groebner_apply!(trace_timing, kat_zp2)

(groebner_time, apply_time)
```

## Batched application over several primes

The paper also sketches the batched `groebner_apply!` workflow. Here the placeholder systems are instantiated explicitly so the example stays executable:

```@example paper_trace
kat_zp3 = map(f -> change_base_ring(GF(2^30 + 9), f), kat)
kat_zp4 = map(f -> change_base_ring(GF(2^30 + 15), f), kat)
kat_zp5 = map(f -> change_base_ring(GF(2^30 + 19), f), kat)

success_batch, (gb2, gb3, gb4, gb5) = groebner_apply!(trace, (kat_zp2, kat_zp3, kat_zp4, kat_zp5))

success_batch && all(isgroebner, (gb2, gb3, gb4, gb5))
```

## Hybrid exact-numeric coefficients

The paper also includes a minimal coefficient type that combines exact modular zero testing with a `BigFloat` payload:

```@example paper_hybrid
import Base: +, -, *, zero, one, inv, iszero, isone

const M = 256
setprecision(BigFloat, M)

struct Hybrid{Zp, Float}
    a::Zp
    b::Float
end

iszero(x::Hybrid) = iszero(x.a)
isone(x::Hybrid) = isone(x.a)
zero(x::Hybrid) = Hybrid(zero(x.a), zero(x.b))
one(x::Hybrid) = Hybrid(one(x.a), one(x.b))

+(x::Hybrid, y::Hybrid) = Hybrid(x.a + y.a, x.b + y.b)
-(x::Hybrid, y::Hybrid) = Hybrid(x.a - y.a, x.b - y.b)
*(x::Hybrid, y::Hybrid) = Hybrid(x.a * y.a, x.b * y.b)
inv(x::Hybrid) = Hybrid(inv(x.a), inv(x.b))

using Groebner, Nemo

to_hybrid(q) = Hybrid(GF(2^30 + 3)(q), BigFloat(q))

sys = Groebner.Examples.hexapod(k = QQ)

sys_exps = map(f -> collect(exponent_vectors(f)), sys)
sys_cfs = map(f -> collect(coefficients(f)), sys)

cfs_hybrid = map(c -> Groebner.CoeffGeneric.(to_hybrid.(c)), sys_cfs)

ring = Groebner.PolyRing(6, DegRevLex(), 0, :generic)
gb_exps, gb_cfs = groebner(ring, sys_exps, cfs_hybrid)

length(gb_exps), length(gb_cfs)
```