## Benchmark results

2024-02-03T15:01:42.740

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.57|
|chandra 12|2.98|
|chandra 13|17.25|
|chandra 14|144.02|
|cyclic 7|0.11|
|cyclic 8|1.26|
|cyclic 9|105.50|
|cyclic 10| - |
|dummy|0.00|
|eco 11|0.37|
|eco 12|2.14|
|eco 13|10.47|
|eco 14|98.45|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|2.24|
|henrion 8|880.36|
|katsura 10|0.77|
|katsura 11|6.83|
|katsura 12|48.55|
|katsura 13|381.16|
|noon 7|0.42|
|noon 8|2.58|
|noon 9|16.02|
|noon 10|145.11|
|noon 11| - |
|reimer 6|0.10|
|reimer 7|0.78|
|reimer 8|23.66|
|reimer 9| - |

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
