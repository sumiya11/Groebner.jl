## Benchmark results

2024-02-09T15:54:14.263

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 12000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|0.47|
|cyclic 8|9.08|
|cyclic 9|5551.77|
|dummy|0.01|
|eco 11|2.15|
|eco 12|12.17|
|eco 13|74.21|
|eco 14|888.17|
|henrion 5|0.03|
|henrion 6|0.35|
|henrion 7|29.50|
|henrion 8| - |
|katsura 10|6.83|
|katsura 11|41.58|
|katsura 12|305.24|
|katsura 13|3367.66|
|noon 7|2.09|
|noon 8|18.18|
|noon 9|201.11|
|reimer 6|0.42|
|reimer 7|6.53|
|reimer 8|257.61|

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
