## Benchmark results

2023-12-29T23:35:20.348

Benchmarked backend: groebner
Benchmark suite: The rationals

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.35|
|cyclic 8|26.73|
|dummy|0.00|
|eco 10|0.64|
|eco 11|3.99|
|eco 12|29.60|
|henrion 6|0.93|
|katsura 10|25.43|
|katsura 11| - |
|katsura 9|4.10|
|noon 8|6.10|
|noon 9|50.79|

*Benchmarking environment:*

* Total RAM (GiB): 188
* Processor: Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
* Julia version: 1.9.2

Versions of the dependencies:

* Primes : 0.5.5
* TimerOutputs : 0.5.23
* PrecompileTools : 1.2.0
* MultivariatePolynomials : 0.5.3
* Combinatorics : 1.0.2
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
