<div align="left">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://github.com/sumiya11/Groebner.jl/raw/master/docs/assets/logo-dark-with-text.svg">
      <img alt="Groebner.jl logo" src="https://github.com/sumiya11/Groebner.jl/raw/master/docs/assets/logo-with-text.svg">
    </picture>
</div>

---

[![Runtests](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sumiya11.github.io/Groebner.jl)
[![codecov](https://codecov.io/github/sumiya11/Groebner.jl/graph/badge.svg?token=J1SZT8ED9S)](https://codecov.io/github/sumiya11/Groebner.jl)

Groebner.jl is a Julia package for computing Groebner bases over fields.
Groebner.jl is distributed under GNU GPL v2.

For documentation please check out https://sumiya11.github.io/Groebner.jl.
For a simple example, see below.

## Installation

You can install Groebner.jl using the Julia package manager. From the Julia REPL, type

```julia
import Pkg; Pkg.add("Groebner")
```

After the installation, you can run package tests with

```julia
import Pkg; Pkg.test("Groebner")
```

## How to use Groebner.jl?

Groebner.jl works with polynomials from AbstractAlgebra.jl, DynamicPolynomials.jl, and Nemo.jl. For example, let's create a ring of polynomials in 3 variables

```julia
using AbstractAlgebra

R, (x1, x2, x3) = QQ["x1", "x2", "x3"]
```

Then, we can define a simple polynomial system

```julia
system = [
  x1 + x2 + x3,
  x1*x2 + x1*x3 + x2*x3,
  x1*x2*x3 - 1
]
```

And compute the Groebner basis by passing the system to `groebner`

```julia
using Groebner

G = groebner(system)
```
```julia
# result
3-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x1 + x2 + x3
 x2^2 + x2*x3 + x3^2
 x3^3 - 1
```

## Contacts

This library is maintained by Alexander Demin (<asdemin_2@edu.hse.ru>).

## Acknowledgement

We would like to acknowledge the developers of the msolve library ([https://msolve.lip6.fr/](https://msolve.lip6.fr/)), as several components of Groebner.jl were adapted from msolve. In our F4 implementation, we adapt and adjust the code of monomial hashtable, critical pair handling and symbolic preprocessing, and linear algebra from msolve. The source code of msolve is available at [https://github.com/algebraic-solving/msolve](https://github.com/algebraic-solving/msolve).

We thank Vladimir Kuznetsov for helpful discussions and providing the sources of his F4 implementation.

We are grateful to The Max Planck Institute for Informatics and The MAX team at l'X for providing computational resources.

## See also

Other software in Julia that can be used for computing Groebner bases:

- [AlgebraicSolving.jl](https://github.com/algebraic-solving/AlgebraicSolving.jl)
- [GroebnerBasis.jl](https://github.com/ederc/GroebnerBasis.jl), deprecated, see [Oscar.jl](https://github.com/oscar-system/Oscar.jl)
- [Singular.jl](https://github.com/oscar-system/Singular.jl)

If you do not see your package here, we either do not know about it, or forgot to include it, sorry! Feel free to open a PR.

## Citing Groebner.jl

If you find Groebner.jl useful in your work, you can star this repository and cite [this paper](https://arxiv.org/abs/2304.06935)

```
@misc{demin2024groebnerjl,
      title={Groebner.jl: A package for Gr\"obner bases computations in Julia}, 
      author={Alexander Demin and Shashi Gowda},
      year={2024},
      eprint={2304.06935},
      archivePrefix={arXiv},
      primaryClass={cs.MS}
}
```
