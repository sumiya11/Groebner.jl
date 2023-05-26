# Groebner.jl

[![Runtests](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sumiya11.github.io/Groebner.jl)


The package provides Groebner bases computation interface in pure Julia with the performance comparable to Singular.
`Groebner.jl` works over finite fields and over the rationals, and supports various monomial orderings.

For documentation and more please check out https://sumiya11.github.io/Groebner.jl

Groebner.jl can be extended : how the code is structured .

## How to use Groebner.jl?

Our package works with polynomials from `AbstractAlgebra.jl`, `DynamicPolynomials.jl`, and `Nemo.jl`. We will demonstrate the usage on a simple example. Lets first create a ring of polynomials in 3 variables

```julia
using AbstractAlgebra
R, (x1, x2, x3) = PolynomialRing(QQ, ["x1", "x2", "x3"]);
```

Then we can define a simple polynomial system

```julia
polys = [
  x1 + x2 + x3,
  x1*x2 + x1*x3 + x2*x3,
  x1*x2*x3 - 1
];
```

And compute the Groebner basis passing the system to `groebner`


```julia
using Groebner
G = groebner(polys)
3-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x1 + x2 + x3
 x2^2 + x2*x3 + x3^2
 x3^3 - 1
```

## Performance

We compare the runtime of our implementation against the ones from `Singular` and `Maple` computer algebra systems. The table below lists measured runtimes of Groebner basis routine for several standard benchmark systems in seconds

|   System    |  Groebner.jl    | Singular | Maple |
| :---:       | :---: | :----: |  :---:   |
| [cyclic-7](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/cyclic12.txt)   | **0.08 s**  | 1.4 s    | **0.08 s** |
| [cyclic-8](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/cyclic13.txt)   |  1.3 s  | 40 s    | **1.1 s** |
| [katsura-10](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/katsura9.txt)    | **0.8 s**  | 71 s    | 0.9 s |
| [katsura-11](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/katsura10.txt)  |  **5.8 s**  | 774 s   | 10 s |
| [eco-12](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/eco10.txt)   |  2.0 s  | 334 s   | **1.6 s** |
| [eco-13](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/eco11.txt)   | **8.8 s**  | 5115 s   | 13 s |
| [noon-7](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/noon7.txt)      |  **0.1 s**  | 0.3 s    | 0.15 s |
| [noon-8](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/systems/standard/noon8.txt)      |  **1.0 s**  | 3.3 s    | 1.1 s |

The bases are computed in `degrevlex` monomial ordering over finite field of characteristic $2^{31}-1$ with all operations single-threaded.

We emphasize that `Groebner.jl` is a specialized library while `Singular` is an extensive general purpose computer algebra system.

## Contacts

This library is maintained by Alexander Demin (<asdemin_2@edu.hse.ru>)

## Acknowledgement

We would like to acknowledge Jérémy Berthomieu, Christian Eder and Mohab Safey El Din as this Library is inspired by their work ["msolve: A Library for Solving Polynomial Systems"](https://arxiv.org/abs/2104.03572). We are also grateful to Max-Planck-Institut für Informatik for assistance in producing benchmarks.

Special thanks goes to Vladimir Kuznetsov for providing the sources of his F4 implementation.

## Citing Groebner.jl

If you find Groebner.jl useful in your work, you can cite [this paper](https://arxiv.org/abs/2304.06935)

```
@misc{groebnerjl2023,
  title = {Groebner.jl: A package for Gr\"obner bases computations in Julia}, 
  author = {Alexander Demin and Shashi Gowda},
  year = {2023},
  eprint = {2304.06935},
  url = {https://arxiv.org/abs/2304.06935}
}
```
