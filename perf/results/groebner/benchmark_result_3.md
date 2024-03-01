## Benchmark results

2024-02-06T11:44:07.983

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.03|
|chandra 10|3.98|
|chandra 11|15.26|
|chandra 12|95.65|
|chandra 13|527.64|
|cyclic 7|0.95|
|cyclic 8|19.67|
|dummy|0.00|
|eco 10|0.48|
|eco 11|2.69|
|eco 12|24.53|
|eco 13|288.08|
|henrion 6|0.47|
|henrion 7|347.42|
|ipp|34.76|
|katsura 9|3.17|
|katsura 10|18.65|
|katsura 11|306.59|
|noon 7|0.81|
|noon 8|4.90|
|noon 9|30.11|
|reimer 6|0.47|
|reimer 7|9.08|
|reimer 8|838.21|

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
