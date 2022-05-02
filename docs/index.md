@def title = "Groebner.jl — Fast Groebner basis computation in Julia"
@def hasmath = false
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->


# Groebner.jl


`Groebner.jl` is a package for fast and generic Gröbner bases computations
based on Faugère's F4 algorithm [[1]](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf) written in pure Julia.

## Installation

To install Groebner.jl, use the Julia package manager:

```julia:install
using Pkg
Pkg.add("Groebner")
```

For a detailed description of exported functions and their input parameters please see page **Interface**. Meanwhile, below are simple usage examples.

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

and compute a *nonreduced* Groebner basis for it
```julia:aagb
using Groebner
basis = groebner(polys, reduced=false)
```

We can also check that a set of polynomials forms a basis
```julia:aaisgb
isgroebner(basis)
```

Now, having a `basis` of the ideal `<polys> = <x^3 + y^2, x*y + x^2>` computed, we can tackle the *ideal membership problem* for `polys`. Recall that a polynomial lies in the ideal iff it's normal form w.r.t a `basis` is zero. Let's check if `x^2y^2 + 2x^3y - xy^2` is a member of `<polys>`!
```julia:aagb
normalform(basis, x^2*y^2 + 2x^3*y - x*y^2)
```

### with `DynamicPolynomials`

We will compute the unique basis for the `noon-2` system, arising from cellular network dynamics analysis [[2]](https://www.jstor.org/stable/2101937):

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
