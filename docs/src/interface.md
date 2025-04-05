# Interface

```@meta
CurrentModule = Groebner
```

## Exported functions

```@docs
groebner
isgroebner
normalform
leading_ideal
dimension
quotient_basis
groebner_with_change_matrix
```

## Monomial orderings

!!! note

    Some frontends, for example, AbstractAlgebra.jl, may not support weighted/product/matrix orderings from Groebner.jl. In such cases, the basis is computed in the ordering requested by user, but the terms of polynomials in the output are ordered w.r.t. some other ordering that is supported by the frontend.

```@docs
Lex
DegLex
DegRevLex
InputOrdering
WeightedOrdering
ProductOrdering
MatrixOrdering
```

## Learn and Apply

```@docs
groebner_learn
groebner_apply!
```
