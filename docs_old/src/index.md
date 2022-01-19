
# Groebner

```@meta
CurrentModule = Groebner
```

`Groebner.jl` is a package for fast and generic Gröbner bases computations
based on Faugère's F4 algorithm [[1]](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf).

## Installation

To install StructuralIdentifiability.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("StructuralIdentifiability")
```

## Interface

```@docs
groebner
```

```@docs
isgroebner
```

```@docs
reduce
```
