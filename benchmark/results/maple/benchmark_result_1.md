## Benchmark results

2024-02-03T15:26:28.713

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.76|
|chandra 12|3.31|
|chandra 13|16.10|
|chandra 14|108.67|
|cyclic 7|0.10|
|cyclic 8|1.35|
|cyclic 9|168.47|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.39|
|eco 12|2.66|
|eco 13|15.26|
|eco 14|122.54|
|henrion 5|0.01|
|henrion 6|0.05|
|henrion 7|3.52|
|henrion 8|1298.38|
|katsura 10|0.95|
|katsura 11|10.30|
|katsura 12|79.35|
|katsura 13|691.95|
|noon 7|0.23|
|noon 8|1.85|
|noon 9|20.14|
|noon 10|198.92|
|noon 11| - |
|reimer 6|0.06|
|reimer 7|1.42|
|reimer 8|43.09|
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
