## Benchmark results

2023-10-03T19:59:45.568

Benchmarked backend: openf4
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 20
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
|henrion 5| - |
|henrion 6| - |
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7| - |
|noon 8| - |
|noon 9| - |
|reimer 6| - |
|reimer 7| - |
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
