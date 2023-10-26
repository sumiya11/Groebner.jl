@def title = "Groebner.jl — Fast Groebner bases in Julia"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

@@im-50
![](/assets/logo-with-text.svg)
@@

Groebner.jl is a package for fast and generic Gröbner bases computations based on Faugère's F4 algorithm[^1] written in Julia.

## Installation

To install Groebner.jl, run the following in the Julia REPL:

```julia
using Pkg; Pkg.add("Groebner")
```

See [Interface](interface) for a description of all exported functions. For a quick introduction to Groebner bases we refer to [Tutorials](tutorial). Meanwhile, below are simple examples.

## Examples

Currently, polynomials from AbstractAlgebra.jl, DynamicPolynomials.jl, and Nemo.jl
are supported as input.

### with AbstractAlgebra.jl

First, import AbstractAlgebra.jl. 
Then, we can create an array of polynomials over a finite field

```julia:install_aa
using Pkg # hide
Pkg.add("AbstractAlgebra") # hide
```

```julia:aaimport
using AbstractAlgebra

R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"])
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

We will compute the basis of the `noon-2` system[^3]

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

We would like to acknowledge Jérémy Berthomieu, Christian Eder, and Mohab Safey El Din as this library is inspired by their work "msolve: A Library for Solving Polynomial Systems"[^2]. We are also grateful to The Max Planck Institute for Informatics and The MAX team at l'X for providing computational resources.

Special thanks goes to Vladimir Kuznetsov for providing the sources of his F4 implementation.

[^1]: [https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf)
[^3]: [https://www.jstor.org/stable/2101937](https://www.jstor.org/stable/2101937)
[^2]: [https://msolve.lip6.fr/](https://msolve.lip6.fr/)
