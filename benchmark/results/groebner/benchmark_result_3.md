## Benchmark results

2024-02-03T00:01:33.383

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 10
- Timeout: 6000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|0.85|
|chandra 10|2.91|
|chandra 11|10.70|
|chandra 12|90.25|
|chandra 13|405.56|
|cyclic 7|0.82|
|cyclic 8|18.80|
|dummy|0.00|
|eco 10|0.59|
|eco 11|2.80|
|eco 12|22.65|
|eco 13|189.21|
|henrion 6|0.49|
|henrion 7|301.23|
|ipp|24.59|
|katsura 9|3.26|
|katsura 10|18.58|
|katsura 11|218.11|
|noon 7|0.61|
|noon 8|4.14|
|noon 9|26.65|
|reimer 6|0.38|
|reimer 7|7.11|
|reimer 8|576.38|

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
