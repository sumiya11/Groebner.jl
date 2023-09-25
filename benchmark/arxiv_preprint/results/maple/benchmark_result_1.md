## Benchmark results

2023-09-25T01:46:00.613

- Benchmarked backend: `maple`
- Workers: 16
- Timeout: 100 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.76|
|cyclic 8|31.19|
|cyclic 9| - |
|dummy|0.01|
|eco 11|39.77|
|eco 12| - |
|eco 13| - |
|henrion 5|0.04|
|henrion 6|2.58|
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|noon 7|13.50|
|noon 8| - |
|noon 9| - |
|reimer 6|0.92|
|reimer 7| - |
|reimer 8| - |

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
* AbstractAlgebra : 0.32.1
* SnoopPrecompile : 1.0.3
