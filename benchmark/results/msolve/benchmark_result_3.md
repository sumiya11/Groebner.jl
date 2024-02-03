## Benchmark results

2024-02-03T03:09:25.084

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.63|
|chandra 10|5.63|
|chandra 11|19.57|
|chandra 12|162.92|
|chandra 13|504.55|
|cyclic 7|1.30|
|cyclic 8|29.24|
|dummy|0.02|
|eco 10|0.60|
|eco 11|4.45|
|eco 12|32.94|
|eco 13|216.88|
|henrion 6|0.76|
|henrion 7|350.72|
|ipp|16.21|
|katsura 9|2.77|
|katsura 10|18.56|
|katsura 11|150.27|
|noon 7|1.38|
|noon 8|11.68|
|noon 9|94.00|
|reimer 6|0.31|
|reimer 7|7.58|
|reimer 8|227.58|

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
