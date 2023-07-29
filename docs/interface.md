@def title = "Groebner.jl — Interface"
@def hasmath = false
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

{{doc isgroebner isgroebner fn}}

{{doc normalform normalform fn}}

<!-- {{doc fglm fglm fn}} -->

{{doc kbase kbase fn}}

## Monomial orderings

A list of all supported monomial orderings.
An ordering can be set by passing the keyword argument `ordering`.
See below for some examples.

*Note that some frontends, for example, `AbstractAlgebra.jl`, may not support weighted/product/matrix orderings. In such cases, the polynomial terms in the output are ordered w.r.t. some other supported ordering.*

{{doc Lex}}

{{doc DegLex}}

{{doc DegRevLex DegRevLex st}}

{{doc InputOrdering InputOrdering st}}

{{doc WeightedOrdering WeightedOrdering st}}

{{doc ProductOrdering ProductOrdering str}}

{{doc MatrixOrdering MatrixOrdering st}}
