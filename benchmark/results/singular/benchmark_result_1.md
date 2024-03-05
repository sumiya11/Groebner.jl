## Benchmark results

2024-02-09T14:05:30.569

Benchmarked backend: singular

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 12000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|1.98|
|cyclic 8|94.22|
|cyclic 9| - |
|dummy|0.00|
|eco 11|58.62|
|eco 12|580.91|
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.40|
|henrion 7|149.01|
|henrion 8| - |
|katsura 10|176.48|
|katsura 11|1170.55|
|katsura 12| - |
|katsura 13| - |
|noon 7|0.47|
|noon 8|4.79|
|noon 9|58.67|
|reimer 6|19.47|
|reimer 7|2633.02|
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
