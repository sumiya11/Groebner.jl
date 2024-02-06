## Benchmark results

2024-02-03T01:01:10.467

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 10
- Timeout: 6500 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|5.60|
|chandra 10|22.04|
|chandra 11|108.12|
|chandra 12|872.18|
|chandra 13|4408.74|
|cyclic 7|1.73|
|cyclic 8|21.95|
|dummy|0.14|
|eco 10|0.78|
|eco 11|4.87|
|eco 12|34.40|
|eco 13|314.34|
|henrion 6|1.67|
|henrion 7|2533.45|
|ipp|139.11|
|katsura 9|8.59|
|katsura 10|78.59|
|katsura 11|1063.73|
|noon 7|1.19|
|noon 8|9.65|
|noon 9|71.94|
|reimer 6|0.70|
|reimer 7|26.74|
|reimer 8|2729.68|

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
