## Benchmark results

2024-02-03T15:43:04.260

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 11|0.57|0.76|0.90|
|chandra 12|2.98|3.31|4.39|
|chandra 13|17.25|16.10|22.34|
|chandra 14|144.02|108.67|185.78|
|cyclic 7|0.11|0.10|0.15|
|cyclic 8|1.26|1.35|1.38|
|cyclic 9|105.50|168.47|125.09|
|cyclic 10| - | - | - |
|dummy|0.00|0.01|0.01|
|eco 11|0.37|0.39|0.54|
|eco 12|2.14|2.66|2.98|
|eco 13|10.47|15.26|18.39|
|eco 14|98.45|122.54|147.85|
|henrion 5|0.00|0.01|0.02|
|henrion 6|0.03|0.05|0.08|
|henrion 7|2.24|3.52|3.52|
|henrion 8|880.36|1298.38|907.76|
|katsura 10|0.77|0.95|1.18|
|katsura 11|6.83|10.30|7.47|
|katsura 12|48.55|79.35|70.53|
|katsura 13|381.16|691.95|483.17|
|noon 7|0.42|0.23|0.23|
|noon 8|2.58|1.85|2.35|
|noon 9|16.02|20.14|20.09|
|noon 10|145.11|198.92|226.64|
|noon 11| - | - | - |
|reimer 6|0.10|0.06|0.07|
|reimer 7|0.78|1.42|1.30|
|reimer 8|23.66|43.09|42.89|
|reimer 9| - | - | - |

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
