## Benchmark results

2024-02-04T04:25:36.935

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.73|
|chandra 12|3.35|
|chandra 13|17.89|
|cyclic 7|0.11|
|cyclic 8|1.41|
|cyclic 9|270.81|
|dummy|0.01|
|eco 11|0.36|
|eco 12|2.23|
|eco 13|14.31|
|eco 14|140.89|
|henrion 5|0.01|
|henrion 6|0.05|
|henrion 7|3.92|
|henrion 8|1190.94|
|katsura 10|1.19|
|katsura 11|10.79|
|katsura 12|69.45|
|katsura 13|922.76|
|katsura 14| - |
|noon 7|0.19|
|noon 8|1.97|
|noon 9|19.40|
|noon 10|209.29|
|reimer 6|0.06|
|reimer 7|1.47|
|reimer 8|47.17|

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
