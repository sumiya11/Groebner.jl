## Benchmark results

2024-01-10T05:33:56.948

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.53|
|cyclic 8|13.38|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.02|
|eco 11|3.60|
|eco 12|21.04|
|eco 13|136.21|
|eco 14| - |
|henrion 5|0.03|
|henrion 6|0.69|
|henrion 7|54.26|
|katsura 10|15.26|
|katsura 11|59.19|
|katsura 12|469.32|
|katsura 13| - |
|noon 7|3.34|
|noon 8|41.71|
|noon 9|482.99|
|noon 10| - |
|reimer 6|0.47|
|reimer 7|13.92|
|reimer 8|501.20|
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
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
