@def title = "Groebner.jl — Interface"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

# Interface

## Exported functions

```julia:load_groebner
using Groebner # hide
```

{{doc groebner groebner fn}}

{{doc groebner_with_change_matrix groebner_with_change_matrix fn}}

{{doc isgroebner isgroebner fn}}

{{doc normalform normalform fn}}

## Monomial orderings

A list of all monomial orderings supported by Groebner.jl.
An ordering can be set by passing it with the keyword argument `ordering`.
See below for some examples.

\note{Some frontends, for example, AbstractAlgebra.jl, may not support weighted/product/matrix orderings from Groebner.jl. In such cases, the basis is computed in the ordering requested by user, but the terms of polynomials in the output are ordered w.r.t. some other ordering that is supported by the frontend.}

{{doc Lex}}

{{doc DegLex}}

{{doc DegRevLex DegRevLex st}}

{{doc InputOrdering InputOrdering st}}

{{doc WeightedOrdering WeightedOrdering st}}

{{doc ProductOrdering ProductOrdering str}}

{{doc MatrixOrdering MatrixOrdering st}}

## Learn and Apply

```julia:load_groebner
using Groebner # hide
```

{{doc groebner_learn groebner_learn fn}}

{{doc groebner_apply! groebner_apply! fn}}

## Low-level interface

```julia:load_groebner
using Groebner # hide
```

Some functions in the interface have a low-level entry point. Low-level functions accept and output ''raw'' exponent vectors and coefficients. This could be convenient when one does not want to depend on a frontend.

For example,

```julia:lowlevel
using Groebner
# define {x y - 1, x^3 + 7 y^2} modulo 65537 in DRL
ring = Groebner.PolyRing(2, Groebner.DegRevLex(), 65537)
monoms = [ [[1, 1], [0, 0]], [[3, 0], [0, 2]] ]
coeffs = [ [1, -1], [1, 7] ]
# compute a GB
gb_monoms, gb_coeffs = Groebner.groebner(ring, monoms, coeffs)
```

The list of functions that provide a low-level entry point: `groebner`, `normalform`, `isgroebner`, `groebner_learn`, `groebner_apply`.

The low-level functions may be faster than their user-facing analogues since they bypass internal checks and conversions. Low-level functions do not make any specific assumptions, that is, all of these are correctly handled in the input: unsorted monomials, nonnormalized coefficients, duplicate terms, aliasing memory.

## Generic coefficients

```julia:load_groebner
using Groebner # hide
```

```julia:install_nemo
using Pkg # hide
Pkg.add("Nemo") # hide
```

Julia grants us the ability to write generic code. One consequence of that for
Groebner.jl is that it can compute Groebner bases over anything that behaves like a field.

For some ground fields Groebner.jl runs an efficient native implementation:
- integers modulo a prime,
- rationals numbers.

For other ground fields, it runs a possibly slower generic fallback. In this case, coefficients of polynomials are treated as black-boxes, which implement field operations: `zero`, `one`, `inv`, `==`, `+`, `*`, `-`.

For example, we can compute a Groebner basis over a univariate rational function field over a finite field:

```julia:generic1
using Groebner, Nemo

R, t = GF(101, 7)["t"]
ff = fraction_field(R)
_, (x, y) = ff["x","y"]

sys = [(t//t+1)*x*y - t^3, y^2 + t]

gb = groebner(sys)
```

Some other functions in Groebner.jl reuse the core F4 algorithm, so they can also be used:

```julia::generic2
@assert isgroebner(gb)
normalform(gb, x*y)
```

### Computing over floating point intervals

Low-level interface supports generic coefficients.

In the following example, we compute a Groebner basis of the `hexapod` system over tuples (Z_p, Interval): each coefficient is treated as a pair, the first coordinate is a finite field element that is used for zero testing, and the second coordinate is a floating point interval with some fixed precision, the payload.

```julia:generic3
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

@info "
    Computing Hexapod over QQ"
c_zp = Groebner.Examples.hexapod(k=AbstractAlgebra.GF(2^30+3));
c_qq = Groebner.Examples.hexapod(k=AbstractAlgebra.QQ);
@time gb_truth = groebner(c_qq);
gbcoeffs_truth = map(f -> collect(coefficients(f)), gb_truth);
@info "
    Coefficient size (in bits): $(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(denominator(c)), f), gbcoeffs_truth))"

@info "
    Computing Hexapod over (Zp, Interval). Precision = $PRECISION bits"
ring = Groebner.PolyRing(nvars(parent(c_qq[1])), Groebner.DegRevLex(), 0);
ring.ground = :generic;

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

Note that if we lower the precision to 256 bits some of the intervals become NaNs.
