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
