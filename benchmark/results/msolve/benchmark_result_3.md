## Benchmark results

2024-01-12T13:17:51.416

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|0.03|
|chandra 10|0.03|
|chandra 11|0.02|
|chandra 12|0.03|
|chandra 13|0.03|
|cyclic 7|1.07|
|cyclic 8|26.84|
|cyclic 9| - |
|dummy|0.03|
|eco 10|0.57|
|eco 11|3.67|
|eco 12|25.40|
|eco 13|443.58|
|henrion 6|3.80|
|henrion 7| - |
|ipp|0.06|
|katsura 9|2.77|
|katsura 10|20.62|
|katsura 11|244.90|
|noon 7|1.56|
|noon 8|12.76|
|noon 9|102.05|
|reimer 6|0.28|
|reimer 7|7.06|
|reimer 8|364.48|

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
