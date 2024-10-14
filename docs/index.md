@def title = "Groebner.jl — Fast Groebner bases in Julia"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

@@im-50
![](/assets/logo-with-text.svg)
@@

Groebner.jl is a package for Gröbner bases computations based on the Faugère's F4 algorithm written in Julia.

## Installation

To install Groebner.jl, run the following in the Julia REPL:

```julia:install-nightly-groebner
using Pkg; # hide
Pkg.add(url="https://github.com/sumiya11/Groebner.jl") # hide
```

```julia
using Pkg; Pkg.add("Groebner")
```

See [Interface](interface) for a description of all exported functions. For a quick introduction to Groebner bases we refer to [Tutorials](tutorial). Meanwhile, below are simple examples.

## Features

Groebner.jl features:

- Gröbner basis computation over integers modulo a prime and over the rationals
- Gröbner trace algorithms
- Multi-threading by default

See [Interface](interface) page for a list of all exported functions.

## Examples

As input, Groebner.jl supports polynomials from AbstractAlgebra.jl, DynamicPolynomials.jl, and Nemo.jl.

### with AbstractAlgebra.jl

First, import AbstractAlgebra.jl. 
Then, we can create an array of polynomials over a finite field

```julia:install_aa
using Pkg # hide
Pkg.add("AbstractAlgebra") # hide
```

```julia:aaimport
using AbstractAlgebra

R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])
polys = [x^2 + y + z, x*y + z];
```

and compute the Groebner basis with the `groebner` command

```julia:aagb
using Groebner

basis = groebner(polys)
```

We can check if a set of polynomials forms a basis

```julia:aaisgb
isgroebner(basis)
```

Groebner.jl also provides several monomial orderings. 
For example, we can eliminate `z` from the above system:

```julia:aagb2
ordering = Lex(z) * DegRevLex(x, y)  # z > x, y
groebner(polys, ordering=ordering)
```

You can find more information on monomial orderings in Groebner.jl in [Monomial Orderings](interface/#monomial_orderings).

### with DynamicPolynomials.jl

We will compute the basis of the `noon-2` system

```julia:install_dynamic
import Pkg; # hide
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

This library is maintained by Alexander Demin ([asdemin_2@edu.hse.ru](mailto:asdemin_2@edu.hse.ru)).

## Acknowledgement

We would like to acknowledge the developers of the msolve library ([https://msolve.lip6.fr/](https://msolve.lip6.fr/)), as several components of Groebner.jl were adapted from msolve. In our F4 implementation, we adapt and adjust the code of monomial hashtable, critical pair handling and symbolic preprocessing, and linear algebra from msolve. The source code of msolve is available at [https://github.com/algebraic-solving/msolve](https://github.com/algebraic-solving/msolve).

We thank Vladimir Kuznetsov for helpful discussions and providing the sources of his F4 implementation.

We are grateful to The Max Planck Institute for Informatics and The MAX team at l'X for providing computational resources.
