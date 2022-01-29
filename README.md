# Groebner.jl

[![Runtests](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/sumiya11/Groebner.jl/actions/workflows/Runtests.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sumiya11.github.io/Groebner.jl)


The package provides Groebner bases computation interface in pure Julia with the performance comparable to Singular.

For documentation and more please check out https://sumiya11.github.io/Groebner.jl

## How to use Groebner.jl?

Our package works with polynomials from `AbstractAlgebra.jl`, `DynamicPolynomials.jl`, and `Nemo.jl`. We will demonstrate the usage on a simple example. Lets first create a ring of polynomials in 3 variables over rationals

```julia
julia> using AbstractAlgebra
julia> R, (x1, x2, x3) = PolynomialRing(QQ, ["x1", "x2", "x3"], ordering=:degrevlex);
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
julia> using Groebner
julia> G = groebner(polys)
3-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x1 + x2 + x3
 x2^2 + x2*x3 + x3^2
 x3^3 - 1
```

## Performance

We compare the runtime of our implementation against the ones from `Singular` and `Maple` computer algebra systems. The table below lists measured runtimes of Groebner basis routine for several standard benchmark systems in seconds

|   System    |  Groebner.jl    | Singular | Maple |
| :---:       | :---: | :----: |  :---:   |
| [cyclic-12](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/cyclic12.txt)   | 0.29s  | **0.01s**    | 0.56s |
| [cyclic-13](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/cyclic13.txt)   |  1.18s  | **0.03s**    | 2.43s |
| [katsura-9](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/katsura9.txt)    | **0.34s**  | 1.13s    | 1.43s |
| [katsura-10](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/katsura10.txt)  |  **2.75s**  | 9.92s   | 5.73s |
| [eco-10](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/eco10.txt)   |  **0.43s**  | 3.22s   | 0.75s |
| [eco-11](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/eco11.txt)   | **3.34s**  | 33.33s   | 3.54s |
| [noon-7](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/noon7.txt)      |  **0.23s**  | 0.40s    | 1.19s|
| [noon-8](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/data/noon8.txt)      |  **2.28s**  | 3.58s    | 8.05s |

The bases are computed in `degrevlex` monomial ordering over finite fields with all operations single-threaded. Systems `cyclic-n` are special in a sense that majority of S-polynomials are redundant, and `Singular` succesfully detects it.

We emphasize that `Groebner.jl` is a specialized library while `Singular`
and `Maple` are extensive general purpose computer algebra systems.

If you discover a system where our package shows bad performance, you are very welcome to submit an issue!  

TODO: speed up rational computations

## Acknowledgement

I would like to acknowledge Jérémy Berthomieu, Christian Eder and Mohab Safey El Din as this Library was adapted from their work ["msolve: A Library for Solving Polynomial Systems"](https://arxiv.org/abs/2104.03572). I also thank Max-Planck-Institut für Informatik for assistance in producing benchmarks.

Special thanks goes to Vladimir Kuznetsov for providing the sources of his F4 implementation.
