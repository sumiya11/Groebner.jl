## Benchmark results

2024-01-26T16:58:11.224

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.25|
|chandra 10|4.53|
|chandra 11|19.83|
|chandra 12|105.81|
|chandra 13| - |
|cyclic 7|1.20|
|cyclic 8|27.28|
|cyclic 9| - |
|dummy|0.02|
|eco 10|0.65|
|eco 11|4.02|
|eco 12|28.25|
|eco 13|314.83|
|henrion 6|0.65|
|henrion 7|404.00|
|ipp|14.18|
|katsura 9|2.21|
|katsura 10|17.29|
|katsura 11|206.39|
|noon 7|1.53|
|noon 8|11.80|
|noon 9|100.94|
|reimer 6|0.30|
|reimer 7|6.14|
|reimer 8|358.83|

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
