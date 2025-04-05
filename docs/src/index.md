# Home

Groebner.jl is a package for computing Gröbner bases written in Julia.

## Installation

To install Groebner.jl, run the following in the Julia REPL:

```julia
using Pkg; Pkg.add("Groebner")
```

## Features

Groebner.jl features:

- Gröbner basis over integers modulo a prime and over the rationals
- Gröbner trace algorithms
- Multi-threading

## Contacts

This library is maintained by Alexander Demin ([asdemin_2@edu.hse.ru](mailto:asdemin_2@edu.hse.ru)).

## Citation

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

## Acknowledgement

We would like to acknowledge the developers of the msolve library ([https://msolve.lip6.fr/](https://msolve.lip6.fr/)), as several components of Groebner.jl were adapted from msolve. In our F4 implementation, we adapt and adjust the code of monomial hashtable, critical pair handling and symbolic preprocessing, and linear algebra from msolve. The source code of msolve is available at [https://github.com/algebraic-solving/msolve](https://github.com/algebraic-solving/msolve).

We thank Vladimir Kuznetsov for helpful discussions and providing the sources of his F4 implementation.

We are grateful to The Max Planck Institute for Informatics, The MAX team at l'X, and the OURAGAN team at Inria for providing computational resources.
