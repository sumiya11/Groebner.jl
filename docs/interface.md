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

{{doc fglm fglm fn}}

{{doc kbase kbase fn}}

## Monomial orderings

A list of supported monomial orderings.
An ordering can be set by passing it as the keyword argument `ordering`.
See below for some examples.

*Note that some frontends (e.g. `AbstractAlgebra.jl`) do not support weighted/block/matrix orderings. In such cases, the polynomial terms in the output may be ordered w.r.t. some other ordering.
Still, the output is a correct Groebner basis in the requested ordering.*

{{doc Lex Lex struct}}

{{doc DegLex DegLex struct}}

{{doc DegRevLex DegRevLex struct}}

{{doc InputOrdering InputOrdering struct}}

{{doc WeightedOrdering WeightedOrdering struct}}

{{doc BlockOrdering BlockOrdering struct}}

{{doc MatrixOrdering MatrixOrdering struct}}
