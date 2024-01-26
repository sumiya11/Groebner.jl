## Benchmark results

2024-01-22T22:02:08.214

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|6.97|
|chandra 12|59.87|
|chandra 13| - |
|chandra 14| - |
|cyclic 7|0.45|
|cyclic 8|11.57|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.02|
|eco 11|3.38|
|eco 12|16.38|
|eco 13|116.34|
|eco 14| - |
|henrion 5|0.04|
|henrion 6|0.37|
|henrion 7|34.52|
|henrion 8| - |
|katsura 10|6.69|
|katsura 11|42.68|
|katsura 12|440.49|
|katsura 13| - |
|noon 7|2.26|
|noon 8|18.99|
|noon 9|431.45|
|noon 10| - |
|noon 11| - |
|reimer 6|0.39|
|reimer 7|11.35|
|reimer 8|379.67|
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
