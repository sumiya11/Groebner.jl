## Benchmark results

2023-09-27T02:33:59.404

Benchmarked backend: msolve
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.13|
|cyclic 8|1.02|
|cyclic 9|61.59|
|dummy|0.02|
|eco 11|0.39|
|eco 12|1.75|
|eco 13|8.55|
|eco 14|70.01|
|henrion 5|0.05|
|henrion 6|0.08|
|henrion 7|2.34|
|katsura 10|0.80|
|katsura 11|4.16|
|katsura 12|25.34|
|katsura 13|151.92|
|noon 7|0.21|
|noon 8|1.23|
|noon 9|9.86|
|reimer 6|0.08|
|reimer 7|0.79|
|reimer 8|16.46|
|reimer 9|796.53|

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
