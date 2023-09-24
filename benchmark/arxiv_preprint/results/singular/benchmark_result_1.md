## Benchmark results

2023-09-25T00:31:37.113

- Benchmarked backend: `singular`
- Workers: 16
- Timeout: 100 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.49|
|cyclic 8|41.28|
|cyclic 9| - |
|dummy|0.00|
|eco 11|31.30|
|eco 12| - |
|eco 13| - |
|henrion 5|0.01|
|henrion 6|0.31|
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|noon 7|0.38|
|noon 8|3.44|
|noon 9|35.32|
|reimer 6|10.01|
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
