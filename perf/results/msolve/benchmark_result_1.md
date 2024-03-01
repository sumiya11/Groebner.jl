## Benchmark results

2024-02-06T17:06:10.405

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|0.15|
|cyclic 8|1.44|
|cyclic 9|270.84|
|dummy|0.01|
|eco 11|0.56|
|eco 12|2.72|
|eco 13|19.39|
|eco 14|259.84|
|henrion 5|0.02|
|henrion 6|0.06|
|henrion 7|3.55|
|henrion 8| - |
|katsura 10|1.20|
|katsura 11|8.90|
|katsura 12|76.77|
|katsura 13|784.47|
|noon 7|0.23|
|noon 8|1.87|
|noon 9|22.32|
|reimer 6|0.08|
|reimer 7|1.23|
|reimer 8|41.38|

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
