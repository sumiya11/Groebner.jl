## Benchmark results

2024-01-12T08:36:22.358

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 5400 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|
|:----|---|---|---|---|
|cyclic 7|0.12|0.19|0.03|0.06|
|cyclic 8|1.44|4.01|0.58|0.75|
|cyclic 9|104.46|643.68|73.29|245.08|
|cyclic 10| - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|
|eco 11|0.37|1.16|0.16|0.28|
|eco 12|2.80|12.76|1.19|1.91|
|eco 13|10.69|133.29|4.89|7.66|
|eco 14|124.62|833.97|66.71|100.08|
|henrion 5|0.00|0.00|0.00|0.00|
|henrion 6|0.07|0.07|0.02|0.05|
|henrion 7|3.08|9.28|1.05|1.98|
|henrion 8| - | - | - | - |
|katsura 10|0.78|5.30|0.31|0.58|
|katsura 11|5.46|74.63|3.36|6.33|
|katsura 12|56.40|592.94|24.49|36.21|
|katsura 13| - | - | - | - |
|noon 7|0.25|0.24|0.04|0.12|
|noon 8|2.05|2.07|0.34|0.74|
|noon 9|15.37|16.85|1.54|2.87|
|noon 10|189.33|194.49|11.82|27.34|
|noon 11| - | - | - | - |
|reimer 6|0.05|0.17|0.01|0.02|
|reimer 7|1.30|1.91|0.18|0.26|
|reimer 8|19.21|58.04|2.55|3.99|
|reimer 9| - | - | - | - |

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
