## Benchmark results

2024-02-07T04:07:14.794

Benchmarked backend: singular

Benchmark suite: The rationals

- Workers: 4
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|2.99|
|chandra 10|44.93|
|chandra 11|341.61|
|chandra 12| - |
|chandra 13| - |
|cyclic 7| - |
|cyclic 8| - |
|dummy|0.00|
|eco 10|38.52|
|eco 11| - |
|eco 12| - |
|eco 13| - |
|henrion 6|6.22|
|henrion 7| - |
|ipp|858.20|
|katsura 9|275.01|
|katsura 10| - |
|katsura 11| - |
|noon 7|1.11|
|noon 8|12.26|
|noon 9|184.11|
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
