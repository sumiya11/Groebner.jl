## Benchmark results

2023-12-31T16:06:12.183

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.12|
|cyclic 8|1.24|
|cyclic 9|137.68|
|cyclic 10| - |
|dummy|0.00|
|eco 11|0.35|
|eco 12|2.08|
|eco 13|9.12|
|eco 14|99.58|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|2.04|
|katsura 10|0.77|
|katsura 11|5.74|
|katsura 12|43.74|
|katsura 13|816.05|
|noon 7|0.19|
|noon 8|1.54|
|noon 9|14.24|
|noon 10|154.35|
|reimer 6|0.05|
|reimer 7|1.23|
|reimer 8|19.93|
|reimer 9| - |

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
