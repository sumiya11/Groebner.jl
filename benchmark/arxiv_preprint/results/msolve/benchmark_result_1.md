## Benchmark results

2023-09-26T03:36:31.168

- Benchmarked backend: `msolve`
- Workers: 8
- Timeout: 12000 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.12|
|cyclic 8|0.98|
|cyclic 9|61.13|
|dummy|0.02|
|eco 11|0.37|
|eco 12|1.70|
|eco 13|8.45|
|eco 14|68.60|
|henrion 5|0.03|
|henrion 6|0.07|
|henrion 7|2.23|
|katsura 10|0.75|
|katsura 11|4.06|
|katsura 12|25.10|
|katsura 13|151.89|
|noon 7|0.19|
|noon 8|1.17|
|noon 9|9.76|
|reimer 6|0.07|
|reimer 7|0.75|
|reimer 8|16.28|
|reimer 9|738.96|

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
