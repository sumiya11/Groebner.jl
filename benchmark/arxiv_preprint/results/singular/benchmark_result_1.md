## Benchmark results

2023-09-27T00:22:25.286

Benchmarked backend: singular
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 100 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|1.45|
|cyclic 8|40.35|
|cyclic 9| - |
|dummy|0.00|
|eco 11|31.38|
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.32|
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.38|
|noon 8|3.44|
|noon 9| - |
|reimer 6|10.12|
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
