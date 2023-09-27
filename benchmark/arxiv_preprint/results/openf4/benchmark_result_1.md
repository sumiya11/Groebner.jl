## Benchmark results

2023-09-27T01:43:48.858

Benchmarked backend: openf4
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 100 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.35|
|cyclic 8|7.31|
|cyclic 9| - |
|dummy|0.01|
|eco 11|1.46|
|eco 12|8.31|
|eco 13| - |
|eco 14| - |
|henrion 5|0.03|
|henrion 6|0.24|
|henrion 7|20.36|
|katsura 10|4.88|
|katsura 11|31.29|
|katsura 12| - |
|katsura 13| - |
|noon 7|1.00|
|noon 8|10.99|
|noon 9| - |
|reimer 6|0.19|
|reimer 7|4.31|
|reimer 8| - |
|reimer 9| - |

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
