# GroebnerBases

Not to forget: we need to change the name

uwu

This repository contains a Julia implementation of the algorithms for Groebner bases computing.

## How to use GroebnerBases.jl?

We will demonstrate the usage on a simple example. Let's first add our package

```julia
julia> import Pkg
Pkg.add(url="https://github.com/sumiya11/GroebnerBases")
```

The package defines several test systems. We will check the algorithm for the set of cyclic polynomials of degree `3` by calling `rootn(3)`

```julia
julia> polys = rootn(3)
3-element Array{AbstractAlgebra.Generic.MPoly{Rational{BigInt}},1}:
 x1 + x2 + x3
 x1*x2 + x1*x3 + x2*x3
 x1*x2*x3 - 1
```

```julia
julia> G = groebner(polys)
julia> G

3-element Array{AbstractAlgebra.Generic.MPoly{Rational{BigInt}},1}:
 x3^3 - 1
 x2^2 + x2*x3 + x3^2
 x1 + x2 + x3
```

Finally, we may want to check if the basis is correct

```julia
julia> is_groebner(G)

true
```
