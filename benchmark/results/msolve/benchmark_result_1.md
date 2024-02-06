## Benchmark results

2024-02-04T04:59:12.529

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.85|
|chandra 12|3.91|
|chandra 13|20.44|
|cyclic 7|0.14|
|cyclic 8|1.38|
|cyclic 9|108.98|
|dummy|0.01|
|eco 11|0.57|
|eco 12|2.57|
|eco 13|12.64|
|eco 14|129.87|
|henrion 5|0.02|
|henrion 6|0.08|
|henrion 7|3.04|
|henrion 8|563.32|
|katsura 10|1.07|
|katsura 11|6.22|
|katsura 12|39.49|
|katsura 13|267.54|
|katsura 14| - |
|noon 7|0.25|
|noon 8|1.73|
|noon 9|16.28|
|noon 10|170.61|
|reimer 6|0.07|
|reimer 7|1.12|
|reimer 8|28.36|

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
