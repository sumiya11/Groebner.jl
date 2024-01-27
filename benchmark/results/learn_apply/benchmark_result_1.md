## Benchmark results

2024-01-27T03:16:11.080

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|Apply 8x, s|
|:----|---|---|---|---|---|
|chandra 11|0.48|1.15|0.05|0.12|0.20|
|chandra 12|2.68|7.36|0.18|0.45|0.69|
|chandra 13|15.06|49.39|0.84|2.16|3.31|
|chandra 14| - | - | - | - | - |
|cyclic 7|0.08|0.18|0.03|0.04|0.06|
|cyclic 8|1.17|3.64|0.50|0.72|1.35|
|cyclic 9| - | - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|0.00|
|eco 11|0.32|1.12|0.14|0.21|0.33|
|eco 12|2.08|9.40|0.85|1.28|2.28|
|eco 13|8.84|88.53|3.58|5.68|12.92|
|eco 14| - | - | - | - | - |
|henrion 5|0.00|0.00|0.00|0.00|0.00|
|henrion 6|0.03|0.05|0.01|0.02|0.06|
|henrion 7|1.93|5.29|0.71|1.13|1.79|
|henrion 8| - | - | - | - | - |
|katsura 10|0.76|4.92|0.30|0.51|0.85|
|katsura 11|4.99|46.82|1.99|3.44|6.24|
|katsura 12| - | - | - | - | - |
|katsura 13| - | - | - | - | - |
|noon 7|0.16|0.19|0.03|0.04|0.08|
|noon 8|1.34|1.56|0.19|0.33|0.46|
|noon 9|12.49|13.06|1.32|2.29|4.06|
|noon 10| - | - | - | - | - |
|reimer 6|0.04|0.08|0.01|0.02|0.02|
|reimer 7|0.77|1.32|0.10|0.16|0.21|
|reimer 8|17.86|42.05|1.70|2.69|4.46|

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
* HostCPUFeatures : 0.1.16
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
