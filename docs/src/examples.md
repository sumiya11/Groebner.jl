# Examples

Groebner.jl supports polynomials from the following frontends:

- AbstractAlgebra.jl
- Nemo.jl
- DynamicPolynomials.jl

Additionally, Groebner.jl has a low-level entry point that accepts raw polynomial data.

## Using AbstractAlgebra.jl

First, we import AbstractAlgebra.jl. 
Then, we create an array of polynomials over a finite field

```@example aa
using AbstractAlgebra

R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])
polys = [x^2 + y + z, x*y + z];
```

and compute a Gröbner basis with the `groebner` command

```@example aa
using Groebner

basis = groebner(polys)
```

We can check if a set of polynomials forms a Gröbner basis

```@example aa
isgroebner(basis)
```

Groebner.jl also provides several monomial orderings. 
For example, we can eliminate `z` from the above system:

```@example aa
ordering = Lex(z) * DegRevLex(x, y)  # z > x, y
groebner(polys, ordering=ordering)
```

You can find more information on monomial orderings in Groebner.jl in [Monomial orderings](@ref).

## Using DynamicPolynomials.jl

Computing the Gröbner basis of some system:

```@example dp
using DynamicPolynomials, Groebner

@polyvar x1 x2
system = [10*x1*x2^2 - 11*x1 + 10,
        10*x1^2*x2 - 11*x2 + 10]

groebner(system)
```

## Using Low-level interface

Some functions in the interface have a low-level entry point. Low-level functions accept and output ''raw'' exponent vectors and coefficients. This could be convenient when one does not want to depend on a frontend.

For example,

```@example
using Groebner

# define {x * y - 1, x^3 + 7 * y^2} modulo 65537 in DRL
ring = Groebner.PolyRing(2, Groebner.DegRevLex(), 65537)
monoms = [ [[1, 1], [0, 0]], [[3, 0], [0, 2]] ]
coeffs = [ [    1,     -1 ], [    1,      7 ] ]

# compute a GB
gb_monoms, gb_coeffs = Groebner.groebner(ring, monoms, coeffs)
```

The list of functions that provide a low-level entry point: `groebner`, `normalform`, `isgroebner`, `groebner_learn`, `groebner_apply`.

Low-level functions may be faster than their user-facing analogues since they bypass data conversions. Low-level functions do not make any specific assumptions on input polynomials, that is, all of these cases are correctly handled: unsorted monomials, non-normalized coefficients, duplicate terms, aliasing memory.

## Generic coefficients

The implementation in Groebner.jl uses a generic type for coefficients. Hence, in theory, Groebner.jl can compute Gröbner bases over any type that behaves like a field.

For the following ground fields Groebner.jl runs an efficient native implementation:
- integers modulo a prime,
- rationals numbers.

For other ground fields, a possibly slower generic fallback is used. In this case, coefficients of polynomials are treated as black-boxes which implement field operations: `zero`, `one`, `inv`, `==`, `+`, `*`, `-`.

For example, we can compute a Gröbner basis over a univariate rational function field over a finite field:

```@example generic2
using Groebner, AbstractAlgebra

R, t = GF(101)["t"]
ff = fraction_field(R)
_, (x, y) = ff["x","y"]

sys = [(t//t+1)*x*y - t^3, y^2 + t]

gb = groebner(sys)
```

Many functions reuse the core implementation, so they can also be used over generic fields:

```julia
@assert isgroebner(gb)
normalform(gb, x*y)
```
### Computing over floating point intervals

In the following example, we combine low-level interface and generic coefficients. 

We are going to compute a basis of the hexapod system over tuples (Z_p, Interval): each coefficient is treated as a pair, the first coordinate is a finite field element used for zero testing, and the second coordinate is a floating point interval with some fixed precision, the payload. For floating point arithmetic, we will be using MPFI.jl.

```@example zp_and_interval
using Pkg;
Pkg.add(url="https://gitlab.inria.fr/ckatsama/mpfi.jl")

import Base: +, -, *, zero, iszero, one, isone, inv
using AbstractAlgebra, Groebner, MPFI

PRECISION = 1024 # For MPFI intervals

struct Zp_And_FloatInterval{Zp, FloatInterval}
    a::Zp
    b::FloatInterval
end

# Pretend it is a field and hakuna matata
+(x::Zp_And_FloatInterval, y::Zp_And_FloatInterval) = Zp_And_FloatInterval(x.a + y.a, x.b + y.b)
*(x::Zp_And_FloatInterval, y::Zp_And_FloatInterval) = Zp_And_FloatInterval(x.a * y.a, x.b * y.b)
-(x::Zp_And_FloatInterval, y::Zp_And_FloatInterval) = Zp_And_FloatInterval(x.a - y.a, x.b - y.b)
zero(x::Zp_And_FloatInterval) = Zp_And_FloatInterval(zero(x.a), zero(x.b))
one(x::Zp_And_FloatInterval) = Zp_And_FloatInterval(one(x.a), one(x.b))
inv(x::Zp_And_FloatInterval) = Zp_And_FloatInterval(inv(x.a), inv(x.b))
iszero(x::Zp_And_FloatInterval) = iszero(x.a)
isone(x::Zp_And_FloatInterval) = isone(x.a)

@info "Computing Hexapod over QQ"
c_zp = Groebner.Examples.hexapod(k=AbstractAlgebra.GF(2^30+3));
c_qq = Groebner.Examples.hexapod(k=AbstractAlgebra.QQ);
@time gb_truth = groebner(c_qq);
gbcoeffs_truth = map(f -> collect(coefficients(f)), gb_truth);
@info "Coefficient size (in bits): $(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(denominator(c)), f), gbcoeffs_truth))"

@info "Computing Hexapod over (Zp, Interval). Precision = $PRECISION bits"
ring = Groebner.PolyRing(nvars(parent(c_qq[1])), Groebner.DegRevLex(), 0, :generic); # Note :generic
exps = map(f -> collect(exponent_vectors(f)), c_zp);
cfs_qq = map(f -> collect(coefficients(f)), c_qq);
cfs_zp = map(f -> collect(coefficients(f)), c_zp);
cfs = map(f -> map(c -> Groebner.CoeffGeneric(Zp_And_FloatInterval(c[1], BigInterval(c[2], precision=PRECISION))), zip(f...)), zip(cfs_zp, cfs_qq));
@time gbexps, gbcoeffs = groebner(ring, exps, cfs);

to_inspect = gbcoeffs[end][end]
@info "
    Inspect one coefficient in the basis:
    Zp         = $(to_inspect.data.a)
    Interval   = $(to_inspect.data.b)
    Diam       = $(diam(to_inspect.data.b))
    Diam (rel) = $(diam_rel(to_inspect.data.b))"

# Sanity check
all_are_inside(x::Zp_And_FloatInterval, truth) = is_inside(BigInterval(truth; precision=PRECISION), x.b)
all_are_inside(x::Groebner.CoeffGeneric, truth) = all_are_inside(x.data, truth)
all_are_inside(x::AbstractVector, truth) = all(map(all_are_inside, x, truth))
@assert all_are_inside(gbcoeffs, gbcoeffs_truth)

# Max |midpoint - truth|
max_error(x::Zp_And_FloatInterval, y; rel=false) = abs(mid(x.b) - y) / ifelse(rel, max(abs(y), 0), 1)
max_error(x::Groebner.CoeffGeneric, y; rel=false) = max_error(x.data, y; rel=rel)
max_error(x::AbstractVector, y::AbstractVector; rel=false) = maximum(map(f -> max_error(f...; rel=rel), zip(x, y)))
@info "
    Max error      : $(max_error(gbcoeffs, gbcoeffs_truth))
    Max error (rel): $(max_error(gbcoeffs, gbcoeffs_truth; rel=true))"

# Max diameter
max_diam(x::Zp_And_FloatInterval; rel=false) = ifelse(rel, diam_rel(x.b), diam(x.b))
max_diam(x::Groebner.CoeffGeneric; rel=false) = max_diam(x.data; rel=rel)
max_diam(x::AbstractVector; rel=false) = maximum(map(f -> max_diam(f; rel=rel), x))
@info "
    Max diam      : $(max_diam(gbcoeffs))
    Max diam (rel): $(max_diam(gbcoeffs; rel=true))"
```

However, if we lower MPFI precision to 256 bits, some of the intervals become NaN.

