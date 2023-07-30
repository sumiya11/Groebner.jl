@def title = "Groebner.jl — Fast Groebner bases in Julia"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->


# Groebner.jl

`Groebner.jl` is a package for fast and generic Gröbner bases computations
based on Faugère's F4 algorithm [[1]](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf) written in pure Julia.

## Installation

To install `Groebner.jl`, execute the following in Julia:

```julia:install
using Pkg; Pkg.add("Groebner")
```

See **Interface** for a description of all exported functions. Find a quick introduction to Groebner bases and a couple of interesting use cases in **Tutorials**. Meanwhile, below are simple examples.

## Examples

Currently, polynomials from `AbstractAlgebra`, `DynamicPolynomials`, and `Nemo`
are supported as input.

### with `AbstractAlgebra`

Let's import `AbstractAlgebra`, create a system over a finite field

```julia:install_aa
Pkg.add("AbstractAlgebra") # hide
```

```julia:aaimport
using AbstractAlgebra

R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"])
polys = [x^3 + y^2, x*y + x^2];
```

and compute the Groebner basis
```julia:aagb
using Groebner

basis = groebner(polys)
```

We can also check that a set of polynomials forms a basis
```julia:aaisgb
isgroebner(basis)
```

We check if $x^2y^2 + 2x^3y - xy^2$ is a member of $\langle x^3 + y^2, xy + x^2 \rangle$, using the `normalform` function

```julia:aagb
normalform(basis, x^2*y^2 + 2x^3*y - x*y^2)
```

### with `DynamicPolynomials`

We will compute the basis of the `noon-2` system [[2]](https://www.jstor.org/stable/2101937):

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

## Contacts

This library is maintained by Alexander Demin (asdemin_2@edu.hse.ru)

## Acknowledgement

We would like to acknowledge Jérémy Berthomieu, Christian Eder and Mohab Safey El Din as this Library is inspired by their work "msolve: A Library for Solving Polynomial Systems". We are also grateful to Max-Planck-Institut für Informatik for assistance in producing benchmarks.

Special thanks goes to Vladimir Kuznetsov for providing the sources of his F4 implementation.
