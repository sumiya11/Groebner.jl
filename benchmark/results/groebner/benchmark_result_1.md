## Benchmark results

2024-02-04T03:52:04.094

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.59|
|chandra 12|3.46|
|chandra 13|18.84|
|cyclic 7|0.11|
|cyclic 8|1.30|
|cyclic 9|150.87|
|dummy|0.00|
|eco 11|0.35|
|eco 12|2.05|
|eco 13|11.75|
|eco 14|133.04|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|2.34|
|henrion 8|871.57|
|katsura 10|0.82|
|katsura 11|6.54|
|katsura 12|50.80|
|katsura 13|483.72|
|katsura 14| - |
|noon 7|0.19|
|noon 8|1.72|
|noon 9|16.14|
|noon 10|188.68|
|reimer 6|0.06|
|reimer 7|1.02|
|reimer 8|25.37|

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
