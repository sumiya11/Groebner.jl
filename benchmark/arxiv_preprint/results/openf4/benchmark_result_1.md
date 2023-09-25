## Benchmark results

2023-09-25T15:48:43.771

- Benchmarked backend: `openf4`
- Workers: 16
- Timeout: 300 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.34|
|cyclic 8|7.62|
|cyclic 9| - |
|dummy|0.01|
|eco 11|1.57|
|eco 12|9.36|
|eco 13|52.36|
|henrion 5|0.02|
|henrion 6|0.23|
|henrion 7|20.73|
|katsura 10|5.07|
|katsura 11|31.50|
|katsura 12|216.15|
|noon 7|0.96|
|noon 8|11.15|
|noon 9|130.39|
|reimer 6|0.18|
|reimer 7|4.23|
|reimer 8|175.34|

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
