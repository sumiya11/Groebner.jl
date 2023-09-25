## Benchmark results

2023-09-25T03:31:55.162

- Benchmarked backend: `groebner`
- Workers: 16
- Timeout: 300 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|4.08|
|cyclic 8|108.93|
|dummy|0.00|
|eco 10|0.78|
|eco 11|6.17|
|henrion 6|3.62|
|katsura 10| - |
|katsura 9|22.24|
|noon 8|8.59|
|noon 9|99.11|

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
