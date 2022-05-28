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
| [cyclic-7](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/cyclic12.txt)   | **0.1s**  | 1.4s    | 0.56s |
| [cyclic-8](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/cyclic13.txt)   |  **2.4s**  | 40s    | 2.43s |
| [katsura-9](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/katsura9.txt)    | **0.10s**  | 1.13s    | 1.43s |
| [katsura-10](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/katsura10.txt)  |  **0.69s**  | 9.92s   | 5.73s |
| [eco-10](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/eco10.txt)   |  **0.17s**  | 3.22s   | 0.75s |
| [eco-11](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/eco11.txt)   | **1.13s**  | 33.33s   | 3.54s |
| [noon-7](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/noon7.txt)      |  **0.20s**  | 0.40s    | 1.19s|
| [noon-8](https://github.com/sumiya11/Groebner.jl/tree/master/benchmark/standard/noon8.txt)      |  **1.88s**  | 3.58s    | 8.05s |

The bases are computed in `degrevlex` monomial ordering over finite field of characteristic $2^{31}-1$ with all operations single-threaded.

We emphasize that `Groebner.jl` is a specialized library while `Singular`
and `Maple` are extensive general purpose computer algebra systems.

If you discover a system where our package shows bad performance, you are very welcome to submit an issue!  

## Acknowledgement

We would like to acknowledge Jérémy Berthomieu, Christian Eder and Mohab Safey El Din as this Library is inspired by their work ["msolve: A Library for Solving Polynomial Systems"](https://arxiv.org/abs/2104.03572). We are also grateful to Max-Planck-Institut für Informatik for assistance in producing benchmarks.

Special thanks goes to Vladimir Kuznetsov for providing the sources of his F4 implementation.
