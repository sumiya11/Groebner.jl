## Benchmark results

2024-01-27T08:58:42.262

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.44|
|chandra 10|5.25|
|chandra 11|21.87|
|chandra 12|115.73|
|chandra 13|709.91|
|cyclic 7|1.11|
|cyclic 8|27.03|
|cyclic 9| - |
|dummy|0.03|
|eco 10|0.66|
|eco 11|4.00|
|eco 12|29.09|
|eco 13|344.04|
|henrion 6|0.63|
|henrion 7|364.54|
|ipp|16.44|
|katsura 9|2.52|
|katsura 10|19.68|
|katsura 11|198.94|
|noon 7|1.75|
|noon 8|11.94|
|noon 9|102.60|
|reimer 6|0.28|
|reimer 7|6.43|
|reimer 8|296.64|

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
