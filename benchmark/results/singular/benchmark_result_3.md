## Benchmark results

2024-01-27T06:50:31.803

Benchmarked backend: singular

Benchmark suite: The rationals

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|3.10|
|chandra 10|35.13|
|chandra 11|284.97|
|chandra 12| - |
|chandra 13| - |
|cyclic 7| - |
|cyclic 8| - |
|cyclic 9| - |
|dummy|0.00|
|eco 10|43.90|
|eco 11| - |
|eco 12| - |
|eco 13| - |
|henrion 6|6.16|
|henrion 7| - |
|ipp|870.29|
|katsura 9|284.81|
|katsura 10| - |
|katsura 11| - |
|noon 7|1.11|
|noon 8|13.23|
|noon 9|172.35|
|reimer 6| - |
|reimer 7| - |
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
* HostCPUFeatures : 0.1.16
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
