## Benchmark results

2024-01-26T16:31:09.203

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|4.84|
|chandra 10|21.17|
|chandra 11|111.70|
|chandra 12|886.62|
|chandra 13| - |
|cyclic 7|1.43|
|cyclic 8|20.50|
|cyclic 9| - |
|dummy|0.14|
|eco 10|0.93|
|eco 11|4.79|
|eco 12|31.59|
|eco 13|352.97|
|henrion 6|2.03|
|henrion 7| - |
|ipp|91.40|
|katsura 9|7.70|
|katsura 10|76.88|
|katsura 11|1218.94|
|noon 7|1.09|
|noon 8|8.15|
|noon 9|66.52|
|reimer 6|0.77|
|reimer 7|24.22|
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
