## Benchmark results

2024-01-27T04:37:58.855

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.28|
|chandra 10|4.20|
|chandra 11|13.27|
|chandra 12|82.00|
|chandra 13|415.78|
|cyclic 7|1.11|
|cyclic 8|27.19|
|cyclic 9| - |
|dummy|0.00|
|eco 10|0.47|
|eco 11|3.16|
|eco 12|26.92|
|eco 13|281.37|
|henrion 6|0.63|
|henrion 7|346.59|
|ipp|93.93|
|katsura 9|3.49|
|katsura 10|24.64|
|katsura 11|269.47|
|noon 7|0.66|
|noon 8|5.37|
|noon 9|36.32|
|reimer 6|0.55|
|reimer 7|10.48|
|reimer 8|661.18|

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
