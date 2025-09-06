## Benchmark results

2023-10-03T20:44:04.983

Benchmarked backend: singular
Benchmark suite: The rationals

- Workers: 20
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7| - |
|cyclic 8| - |
|dummy|0.00|
|eco 10|29.41|
|eco 11| - |
|henrion 6|4.05|
|katsura 10| - |
|katsura 9|127.50|
|noon 8|10.66|
|noon 9|109.81|

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
