## Benchmark results

2023-09-27T01:22:28.991

Benchmarked backend: groebner
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 10 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7| - |
|cyclic 8| - |
|cyclic 9| - |
|dummy| - |
|eco 11| - |
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5|0.00|
|henrion 6|0.04|
|henrion 7|1.73|
|katsura 10|0.69|
|katsura 11|4.48|
|katsura 12|31.93|
|katsura 13|216.65|
|noon 7|0.17|
|noon 8|1.58|
|noon 9|15.71|
|reimer 6|0.04|
|reimer 7|0.78|
|reimer 8|16.92|

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
