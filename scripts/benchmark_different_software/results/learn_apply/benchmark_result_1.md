## Benchmark results

2024-02-09T12:36:24.270

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 10000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|Apply 8x, s|
|:----|---|---|---|---|---|
|cyclic 7|0.11|0.18|0.03|0.06|3.09|
|cyclic 8|1.20|3.75|0.53|0.76|4.48|
|cyclic 9|98.96|691.77|75.67|141.23|277.85|
|dummy|0.00|0.00|0.00|0.00|2.99|
|eco 11|0.33|1.19|0.14|0.22|3.64|
|eco 12|2.05|9.38|0.91|1.40|5.92|
|eco 13|9.26|98.25|4.07|6.49|15.39|
|eco 14|97.60|926.51|55.31|104.95|287.67|
|henrion 5|0.00|0.03|0.00|0.00|3.15|
|henrion 6|0.03|0.05|0.04|0.02|3.27|
|henrion 7|2.13|5.95|0.89|1.56|5.71|
|henrion 8| - | - | - | - | - |
|katsura 10|0.82|5.49|0.32|0.61|4.37|
|katsura 11|5.66|49.95|2.19|3.84|11.35|
|katsura 12|45.52|523.29|15.96|28.73|62.91|
|katsura 13| - | - | - | - | - |
|noon 7|0.21|0.24|0.03|0.07|3.59|
|noon 8|1.66|1.68|0.19|0.39|4.49|
|noon 9|14.31|15.14|1.46|2.81|9.33|
|reimer 6|0.08|0.06|0.01|0.02|3.45|
|reimer 7|0.89|1.48|0.11|0.18|3.93|
|reimer 8|19.65|44.75|1.84|2.97|8.90|

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
