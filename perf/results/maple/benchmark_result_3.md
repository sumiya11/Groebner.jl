## Benchmark results

2024-02-06T12:16:34.827

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|5.37|
|chandra 10|25.70|
|chandra 11|129.58|
|chandra 12|1052.14|
|chandra 13| - |
|cyclic 7|1.41|
|cyclic 8|23.82|
|dummy|0.01|
|eco 10|0.71|
|eco 11|4.87|
|eco 12|35.06|
|eco 13|496.48|
|henrion 6|1.74|
|henrion 7| - |
|ipp|173.43|
|katsura 9|9.00|
|katsura 10|84.81|
|katsura 11|1318.31|
|noon 7|1.03|
|noon 8|9.23|
|noon 9|98.78|
|reimer 6|0.57|
|reimer 7|32.33|
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
