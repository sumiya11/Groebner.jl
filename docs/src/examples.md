# Examples

Groebner.jl supports polynomials from the following frontends:

- AbstractAlgebra.jl
- Nemo.jl
- DynamicPolynomials.jl

Additionally, Groebner.jl provides a low-level entry point that accepts raw polynomial data.

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

```@example generic-2
using Groebner, AbstractAlgebra

R, t = GF(101)["t"]
ff = fraction_field(R)
_, (x, y) = ff["x","y"]

sys = [(t//t+1)*x*y - t^3, y^2 + t]

gb = groebner(sys)
```

Many functions reuse the core implementation, so they can also be used over generic fields:

```@example generic-2
@assert isgroebner(gb)
normalform(gb, x*y)
```
