## Benchmark results

2023-09-26T00:12:55.748

- Benchmarked backend: `groebner`
- Workers: 8
- Timeout: 12000 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.08|
|cyclic 8|1.15|
|cyclic 9|95.93|
|dummy|0.00|
|eco 11|0.30|
|eco 12|1.80|
|eco 13|8.31|
|eco 14|79.48|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|1.68|
|katsura 10|0.65|
|katsura 11|4.49|
|katsura 12|31.71|
|katsura 13|218.26|
|noon 7|0.14|
|noon 8|1.60|
|noon 9|15.73|
|reimer 6|0.04|
|reimer 7|0.63|
|reimer 8|16.53|
|reimer 9|698.89|

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
