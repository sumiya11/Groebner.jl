## Benchmark results

2024-02-06T12:27:10.131

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.29|
|chandra 10|4.55|
|chandra 11|21.31|
|chandra 12|104.74|
|chandra 13|553.96|
|cyclic 7|1.13|
|cyclic 8|26.10|
|dummy|0.02|
|eco 10|0.59|
|eco 11|3.72|
|eco 12|29.19|
|eco 13|269.46|
|henrion 6|0.67|
|henrion 7|391.25|
|ipp|16.87|
|katsura 9|2.27|
|katsura 10|17.54|
|katsura 11|167.57|
|noon 7|1.51|
|noon 8|11.37|
|noon 9|91.26|
|reimer 6|0.32|
|reimer 7|6.68|
|reimer 8|256.74|

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
