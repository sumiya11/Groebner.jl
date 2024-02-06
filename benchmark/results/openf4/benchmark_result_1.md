## Benchmark results

2024-02-04T04:43:53.759

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|14.24|
|chandra 12|53.13|
|chandra 13|274.22|
|cyclic 7|0.48|
|cyclic 8|9.64|
|cyclic 9| - |
|dummy|0.01|
|eco 11|2.26|
|eco 12|12.48|
|eco 13|98.25|
|eco 14|716.15|
|henrion 5|0.04|
|henrion 6|0.63|
|henrion 7|41.39|
|henrion 8| - |
|katsura 10|8.65|
|katsura 11|61.22|
|katsura 12|567.10|
|katsura 13| - |
|katsura 14| - |
|noon 7|1.69|
|noon 8|20.70|
|noon 9|302.98|
|noon 10| - |
|reimer 6|0.32|
|reimer 7|9.40|
|reimer 8|373.76|

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
