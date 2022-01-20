@def title = "Groebner.jl — Fast Groebner basis computation in Julia"
@def hasmath = false
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->


# Groebner.jl


`Groebner.jl` is a package for fast and generic Gröbner bases computations
based on Faugère's F4 algorithm [[1]](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf).

## Installation

To install Groebner.jl, use the Julia package manager:

```julia:install
using Pkg
Pkg.add(url="https://github.com/sumiya11/Groebner.jl")
```

For a detailed description of exported functions and their input parameters please see **Interface**. Meanwhile, below are simple usage examples.

## Examples

Currently, polynomials from `AbstractAlgebra`, `DynamicPolynomials`, and `Nemo`
are supported as input.

### with `AbstractAlgebra`

Let's import `AbstractAlgebra`, create a system over a finite field..

```julia:install_aa
Pkg.add("AbstractAlgebra") # hide
```

```julia:aaimport
using AbstractAlgebra
R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"])
polys = [x^3 + y^2, x*y + x^2];
```

..compute a *nonreduced* Groebner basis for it..
```julia:aagb
basis = groebner(polys, reduced=false)
```

..and check correctness of result
```julia:aaisgb
isgroebner(basis)
```

### with `DynamicPolynomials`

We will compute the basis for the `noon-2` system, arising from cellular network dynamics analysis [[2]](https://www.jstor.org/stable/2101937):

```julia:install_dynamic
Pkg.add("DynamicPolynomials") # hide
```

```julia:aaimport
using DynamicPolynomials
@polyvar x1 x2
system = [10*x1*x2^2 - 11*x1 + 10,
        10*x1^2*x2 - 11*x2 + 10]

groebner(system)
```
