## Benchmark results

2023-09-25T03:34:56.555

- Benchmarked backend: `msolve`
- Workers: 16
- Timeout: 300 s

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|2.03|
|cyclic 8|42.60|
|dummy|0.03|
|eco 10|0.57|
|eco 11|3.58|
|henrion 6|4.28|
|katsura 10|22.80|
|katsura 9|3.03|
|noon 8|10.07|
|noon 9|86.66|

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
