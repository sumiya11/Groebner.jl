## Benchmark results

2023-09-27T07:21:53.204

Benchmarked backend: groebner
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.09|
|cyclic 8|1.25|
|cyclic 9|96.04|
|dummy|0.00|
|eco 11|0.32|
|eco 12|1.81|
|eco 13|8.53|
|eco 14|78.01|
|henrion 5|0.00|
|henrion 6|0.05|
|henrion 7|1.68|
|katsura 10|0.68|
|katsura 11|4.59|
|katsura 12|32.09|
|katsura 13|219.01|
|noon 7|0.16|
|noon 8|1.65|
|noon 9|16.47|
|reimer 6|0.04|
|reimer 7|0.66|
|reimer 8|17.22|
|reimer 9|681.04|

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
