## Benchmark results

2023-09-25T16:17:24.805

- Benchmarked backend: `groebner`
- Workers: 16
- Timeout: 300 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.10|
|cyclic 8|1.18|
|cyclic 9|96.56|
|dummy|0.00|
|eco 11|0.32|
|eco 12|1.84|
|eco 13|8.25|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|1.75|
|katsura 10|0.67|
|katsura 11|4.50|
|katsura 12|31.70|
|noon 7|0.16|
|noon 8|1.52|
|noon 9|15.80|
|reimer 6|0.05|
|reimer 7|0.71|
|reimer 8|16.84|

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
