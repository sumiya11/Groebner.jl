## Benchmark results

2024-02-06T20:39:54.517

Benchmarked backend: singular

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 7200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|1.95|
|cyclic 8|88.30|
|cyclic 9| - |
|dummy|0.00|
|eco 11|57.00|
|eco 12|628.17|
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.41|
|henrion 7|94.75|
|henrion 8| - |
|katsura 10|107.25|
|katsura 11|1387.56|
|katsura 12| - |
|katsura 13| - |
|noon 7|0.49|
|noon 8|4.93|
|noon 9|47.89|
|reimer 6|19.51|
|reimer 7|3106.54|
|reimer 8| - |

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
