# FastGroebner.jl

[![Runtests](https://github.com/sumiya11/FastGroebner.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/FastGroebner.jl/actions/workflows/Runtests.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sumiya11.github.io/FastGroebner.jl/dev)


The package provides Groebner bases computation interface with the *performance
comparable to Singular*.

**Please note that the package is still very WIP.**

## How to use FastGroebner.jl?

We will demonstrate the usage on a simple example. Lets first import `AbstractAlgebra`
and create a ring of polynomials in 3 variables over rationals

```julia
julia> using AbstractAlgebra
julia> R, (x1, x2, x3) = PolynomialRing(GF(2^31 - 1), ["x1", "x2", "x3"], ordering=:degrevlex);
```

Then we can define a simple cyclic polynomial system

```julia
julia> polys = [
  x1 + x2 + x3,
  x1*x2 + x1*x3 + x2*x3,
  x1*x2*x3 - 1
];
```

And compute the Groebner basis passing the system to `groebner`


```julia
julia> using FastGroebner
julia> G = groebner(polys)
3-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x1 + x2 + x3
 x2^2 + x2*x3 + x3^2
 x3^3 - 1
```

## Performance

We compare the runtime of our implementation against the one of Singular Groebner bases routine. The table below lists measured runtimes for several standard benchmark systems in seconds

|   System    | Size  | Ours    | Singular |
| :---:       | :---: |  :----: |  :---:   |
| cyclic-12   |  23   |  0.18s  | **0.03s**    |
| cyclic-13   |  25   |  0.58s  | **0.10s**    |
| katsura-9   |  145   |  **0.97s**  | 2.81s    |
| katsura-10  |  274   |  **7.82s**  | 25.426s    |
| noon-6      |  187   |  **0.09s**  | 0.11s    |
| noon-7      |  495   |  **0.71s**  | 0.89s    |

The bases are computed in `degrevlex` monomial ordering over finite fields having all operations single-threaded. Benchmarks sources are located in the `benchmark` folder.

## Acknowledgement

I would like to acknowledge Jérémy Berthomieu, Christian Eder and Mohab Safey El Din as this Library was adapted from their work ["msolve: A Library for Solving Polynomial Systems"](https://arxiv.org/abs/2104.03572). I also thank Vladimir Kuznetsov for providing the sources of his F4 implementation.
