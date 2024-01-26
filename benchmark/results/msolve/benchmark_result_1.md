## Benchmark results

2024-01-26T15:22:43.818

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.88|
|chandra 12|3.84|
|chandra 13|20.56|
|chandra 14|133.08|
|cyclic 7|0.15|
|cyclic 8|1.36|
|cyclic 9|118.97|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.52|
|eco 12|2.59|
|eco 13|13.42|
|eco 14|127.94|
|henrion 5|0.02|
|henrion 6|0.06|
|henrion 7|3.26|
|henrion 8| - |
|katsura 10|1.08|
|katsura 11|6.15|
|katsura 12|43.52|
|katsura 13| - |
|noon 7|0.22|
|noon 8|1.88|
|noon 9|16.01|
|noon 10| - |
|noon 11| - |
|reimer 6|0.08|
|reimer 7|1.21|
|reimer 8|30.17|
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
