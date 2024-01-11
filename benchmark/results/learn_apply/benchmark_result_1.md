## Benchmark results

2024-01-11T09:17:31.780

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|
|-----|---|---|---|---|
|cyclic 7|0.10|0.24|0.03|0.05|
|cyclic 8|1.37|4.24|0.67|0.93|
|cyclic 9| - | - | - | - |
|cyclic 10| - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|
|eco 11|0.36|1.15|0.16|0.27|
|eco 12|2.85|11.77|1.08|1.73|
|eco 13|14.58|144.95|8.58|13.82|
|eco 14| - | - | - | - |
|henrion 5|0.00|0.00|0.00|0.00|
|henrion 6|0.07|0.06|0.02|0.05|
|henrion 7|2.72|9.00|0.93|1.47|
|henrion 8| - | - | - | - |
|katsura 10|0.98|5.58|0.34|0.70|
|katsura 11|10.15|101.82|3.81|7.26|
|katsura 12| - | - | - | - |
|katsura 13| - | - | - | - |
|noon 7|0.19|0.34|0.04|0.11|
|noon 8|2.02|2.37|0.31|0.62|
|noon 9|18.01|18.41|1.88|3.49|
|noon 10| - | - | - | - |
|reimer 6|0.06|0.07|0.06|0.02|
|reimer 7|1.30|1.71|0.12|0.24|
|reimer 8|29.25|84.44|2.58|4.25|
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
