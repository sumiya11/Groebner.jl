## Benchmark results

2024-01-27T02:30:44.878

Benchmarked backend: groebner

Benchmark suite: Integers modulo 1031

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.48|
|chandra 12|2.43|
|chandra 13|20.02|
|chandra 14| - |
|cyclic 7|0.08|
|cyclic 8|1.19|
|cyclic 9| - |
|dummy|0.00|
|eco 11|0.31|
|eco 12|1.84|
|eco 13|8.61|
|eco 14|163.92|
|henrion 5|0.00|
|henrion 6|0.04|
|henrion 7|2.10|
|henrion 8| - |
|katsura 10|0.77|
|katsura 11|5.15|
|katsura 12|52.61|
|katsura 13| - |
|noon 7|0.18|
|noon 8|1.32|
|noon 9|12.09|
|noon 10| - |
|reimer 6|0.04|
|reimer 7|0.75|
|reimer 8|18.54|

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
