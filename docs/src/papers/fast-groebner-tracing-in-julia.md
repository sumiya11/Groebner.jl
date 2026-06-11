# Groebner.jl: Fast Gröbner Tracing in Julia

This page mirrors the code snippets from the paper _Groebner.jl: Fast Gröbner Tracing in Julia_ that appeared in the proceedings of the [International Congress on Mathematical Software 2026](https://icms-conference.org/2026/).

## Setup

```julia
using Pkg; Pkg.add("Groebner"); Pkg.add("Nemo")
```

## Generic coefficient arithmetic

```@example paper_generic_coefficients
using Groebner, Nemo

K, t = rational_function_field(QQ, "t")
R, (x, y, z) = K["x", "y", "z"]

F = [t * x^2 * y + y^3, x * y^2 + 3 // t * z]

G = groebner(F);
(length(G), isgroebner(G))
```

## Learn and apply

```@example paper_trace
using Groebner, Nemo

kat = Groebner.Examples.katsuran(10)

kat_zp1 = map(f -> change_base_ring(GF(2^30 + 3), f), kat)
kat_zp2 = map(f -> change_base_ring(GF(2^30 + 7), f), kat)

trace, _ = groebner_learn(kat_zp1)
success, gb = groebner_apply!(trace, kat_zp2)

@assert success && isgroebner(gb)
length(gb)
```

## Timing the apply phase

```@example paper_trace
using BenchmarkTools
@btime groebner(kat_zp2);
```

```@example paper_trace
@btime groebner_apply!(trace, kat_zp2);
```

## Batched application over several primes

```@example paper_trace
kat_zp3 = map(f -> change_base_ring(GF(2^30 + 9), f), kat)
kat_zp4 = map(f -> change_base_ring(GF(2^30 + 15), f), kat)
kat_zp5 = map(f -> change_base_ring(GF(2^30 + 19), f), kat)

success, (gb2, gb3, gb4, gb5) = groebner_apply!(trace, (kat_zp2, kat_zp3, kat_zp4, kat_zp5))

@assert success && all(isgroebner, (gb2, gb3, gb4, gb5))
map(length, (gb2, gb3, gb4, gb5))
```

## Hybrid exact-numeric coefficients

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
```

```@example paper_hybrid
using Groebner, Nemo

to_hybrid(q) = Hybrid(GF(2^30 + 3)(q), BigFloat(q))

sys = Groebner.Examples.hexapod(k = QQ)

sys_exps = map(f -> collect(exponent_vectors(f)), sys)
sys_cfs = map(f -> collect(coefficients(f)), sys)

cfs_hybrid = map(c -> Groebner.CoeffGeneric.(to_hybrid.(c)), sys_cfs)

ring = Groebner.PolyRing(6, DegRevLex(), 0, :generic)
gb_exps, gb_cfs = groebner(ring, sys_exps, cfs_hybrid);
```
