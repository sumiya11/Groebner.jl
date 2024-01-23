## Benchmark results

2024-01-23T14:39:51.773

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.52|
|chandra 12|2.76|
|chandra 13|16.85|
|chandra 14|124.44|
|cyclic 7|0.11|
|cyclic 8|1.26|
|cyclic 9|163.44|
|cyclic 10| - |
|dummy|0.00|
|eco 11|0.32|
|eco 12|2.17|
|eco 13|9.18|
|eco 14|97.59|
|henrion 5|0.00|
|henrion 6|0.04|
|henrion 7|2.14|
|henrion 8| - |
|katsura 10|0.80|
|katsura 11|5.86|
|katsura 12|36.37|
|katsura 13| - |
|noon 7|0.16|
|noon 8|1.54|
|noon 9|14.00|
|noon 10|220.48|
|noon 11| - |
|reimer 6|0.05|
|reimer 7|0.95|
|reimer 8|19.32|
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
