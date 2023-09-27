## Benchmark results

2023-09-27T02:24:08.357

Benchmarked backend: groebner
Benchmark suite: The rationals

- Workers: 16
- Timeout: 600 s
- Aggregated over: 2 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|4.15|
|cyclic 8|106.70|
|dummy|0.00|
|eco 10|0.77|
|eco 11|6.11|
|henrion 6|3.52|
|katsura 10| - |
|katsura 9|22.26|
|noon 8|8.67|
|noon 9|99.24|

*Benchmarking environment:*

* Total RAM (GiB): 2003
* Processor: AMD EPYC 7702 64-Core Processor                
* Julia version: 1.9.1

Versions of the dependencies:

* Combinatorics : 1.0.2
* MultivariatePolynomials : 0.5.2
* Primes : 0.5.4
* ExprTools : 0.1.10
* SIMD : 3.4.5
* AbstractAlgebra : 0.32.3
* SnoopPrecompile : 1.0.3
