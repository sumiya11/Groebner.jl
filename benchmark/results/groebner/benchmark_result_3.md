## Benchmark results

2024-01-22T23:54:07.650

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 16
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.22|
|chandra 10|2.92|
|chandra 11|10.60|
|chandra 12|87.26|
|chandra 13| - |
|cyclic 7|0.88|
|cyclic 8|21.57|
|cyclic 9| - |
|dummy|0.00|
|eco 10|0.40|
|eco 11|2.88|
|eco 12|20.95|
|eco 13| - |
|henrion 6|0.53|
|henrion 7| - |
|ipp|86.51|
|katsura 9|2.70|
|katsura 10|20.97|
|katsura 11| - |
|noon 7|0.48|
|noon 8|3.95|
|noon 9|26.13|
|reimer 6|0.42|
|reimer 7|7.46|
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
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
