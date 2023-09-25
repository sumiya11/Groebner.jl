## Benchmark results

2023-09-25T17:57:16.943

- Benchmarked backend: `msolve`
- Workers: 16
- Timeout: 300 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.13|
|cyclic 8|0.99|
|cyclic 9|61.38|
|dummy|0.04|
|eco 11|0.38|
|eco 12|1.72|
|eco 13|8.49|
|henrion 5|0.04|
|henrion 6|0.07|
|henrion 7|2.30|
|katsura 10|0.86|
|katsura 11|4.16|
|katsura 12|25.35|
|noon 7|0.19|
|noon 8|1.28|
|noon 9|9.57|
|reimer 6|0.06|
|reimer 7|0.81|
|reimer 8|16.06|

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
