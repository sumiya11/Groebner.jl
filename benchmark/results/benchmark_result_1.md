## Benchmark results

2024-02-06T17:06:11.759

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|cyclic 7|0.13|0.10|0.15|
|cyclic 8|1.46|1.23|1.44|
|cyclic 9|259.19|340.73|270.84|
|dummy|0.00|0.01|0.01|
|eco 11|0.38|0.35|0.56|
|eco 12|2.23|2.00|2.72|
|eco 13|10.98|13.18|19.39|
|eco 14|213.55|216.66|259.84|
|henrion 5|0.00|0.01|0.02|
|henrion 6|0.09|0.05|0.06|
|henrion 7|2.83|3.15|3.55|
|henrion 8| - | - | - |
|katsura 10|0.94|1.41|1.20|
|katsura 11|9.20|7.72|8.90|
|katsura 12|80.44|52.66|76.77|
|katsura 13| - | - |784.47|
|noon 7|0.21|0.23|0.23|
|noon 8|2.00|1.95|1.87|
|noon 9|19.49|18.67|22.32|
|reimer 6|0.05|0.06|0.08|
|reimer 7|1.07|1.16|1.23|
|reimer 8|31.75|31.81|41.38|

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
