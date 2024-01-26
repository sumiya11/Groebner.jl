## Benchmark results

2024-01-23T15:08:36.354

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 11|0.52|0.62|0.93|
|chandra 12|2.76|2.72|3.90|
|chandra 13|16.85|14.51|23.23|
|chandra 14|124.44|103.12|152.35|
|cyclic 7|0.11|0.10|0.14|
|cyclic 8|1.26|1.20|1.43|
|cyclic 9|163.44|209.38|151.28|
|cyclic 10| - | - | - |
|dummy|0.00|0.01|0.01|
|eco 11|0.32|0.36|0.54|
|eco 12|2.17|1.95|2.81|
|eco 13|9.18|10.94|16.03|
|eco 14|97.59|125.95|150.29|
|henrion 5|0.00|0.01|0.02|
|henrion 6|0.04|0.05|0.08|
|henrion 7|2.14|3.03|3.14|
|henrion 8| - | - | - |
|katsura 10|0.80|1.26|1.20|
|katsura 11|5.86|7.33|9.45|
|katsura 12|36.37|57.21|58.58|
|katsura 13| - |987.88| - |
|noon 7|0.16|0.19|0.23|
|noon 8|1.54|2.12|1.90|
|noon 9|14.00|15.12|21.43|
|noon 10|220.48|238.63|266.79|
|noon 11| - | - | - |
|reimer 6|0.05|0.06|0.07|
|reimer 7|0.95|1.67|1.18|
|reimer 8|19.32|28.65|41.91|
|reimer 9| - | - | - |

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
