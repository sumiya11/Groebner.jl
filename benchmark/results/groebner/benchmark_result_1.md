## Benchmark results

2023-10-13T14:13:57.207

Benchmarked backend: groebner
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 60 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7| - |
|cyclic 8|2.59|
|cyclic 9| - |
|dummy| - |
|eco 11| - |
|eco 12| - |
|eco 13| - |
|henrion 5| - |
|henrion 6|0.12|
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.58|
|noon 8| - |
|noon 9| - |
|reimer 6| - |
|reimer 7|1.88|
|reimer 8| - |
|reimer 9| - |

*Benchmarking environment:*

* Total RAM (GiB): 7
* Processor: Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz
* Julia version: 1.9.1

Versions of the dependencies:

* Combinatorics : 1.0.2
* MultivariatePolynomials : 0.5.2
* Primes : 0.5.4
* ExprTools : 0.1.10
* SIMD : 3.4.5
* AbstractAlgebra : 0.32.3
* SnoopPrecompile : 1.0.3
