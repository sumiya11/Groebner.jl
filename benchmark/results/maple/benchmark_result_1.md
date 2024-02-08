## Benchmark results

2024-02-06T16:48:54.678

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|0.10|
|cyclic 8|1.23|
|cyclic 9|340.73|
|dummy|0.01|
|eco 11|0.35|
|eco 12|2.00|
|eco 13|13.18|
|eco 14|216.66|
|henrion 5|0.01|
|henrion 6|0.05|
|henrion 7|3.15|
|henrion 8| - |
|katsura 10|1.41|
|katsura 11|7.72|
|katsura 12|52.66|
|katsura 13| - |
|noon 7|0.23|
|noon 8|1.95|
|noon 9|18.67|
|reimer 6|0.06|
|reimer 7|1.16|
|reimer 8|31.81|

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
