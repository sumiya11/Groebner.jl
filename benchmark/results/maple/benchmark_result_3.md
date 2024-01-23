## Benchmark results

2024-01-23T00:11:24.086

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 16
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|4.87|
|chandra 10|21.35|
|chandra 11|121.88|
|chandra 12| - |
|chandra 13| - |
|cyclic 7|1.37|
|cyclic 8|21.14|
|cyclic 9| - |
|dummy|0.01|
|eco 10|0.66|
|eco 11|4.64|
|eco 12|31.94|
|eco 13|520.63|
|henrion 6|1.69|
|henrion 7| - |
|ipp|148.38|
|katsura 9|7.76|
|katsura 10|76.62|
|katsura 11| - |
|noon 7|1.03|
|noon 8|8.21|
|noon 9|67.71|
|reimer 6|0.65|
|reimer 7|25.56|
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
