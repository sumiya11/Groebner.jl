## Benchmark results

2024-01-23T01:56:01.563

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 16
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.50|
|chandra 10|4.51|
|chandra 11|21.99|
|chandra 12|118.94|
|chandra 13| - |
|cyclic 7|1.37|
|cyclic 8|26.78|
|cyclic 9| - |
|dummy|0.03|
|eco 10|0.62|
|eco 11|3.85|
|eco 12|28.75|
|eco 13|380.05|
|henrion 6|0.63|
|henrion 7|428.35|
|ipp|20.31|
|katsura 9|2.77|
|katsura 10|18.49|
|katsura 11|223.14|
|noon 7|1.60|
|noon 8|12.65|
|noon 9|110.47|
|reimer 6|0.31|
|reimer 7|6.69|
|reimer 8|339.74|

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
