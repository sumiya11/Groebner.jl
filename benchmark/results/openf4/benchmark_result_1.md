## Benchmark results

2023-12-30T23:01:44.844

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.46|
|cyclic 8|9.53|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.01|
|eco 11|2.36|
|eco 12|18.27|
|eco 13|76.21|
|eco 14| - |
|henrion 5|0.04|
|henrion 6|0.33|
|henrion 7|28.49|
|katsura 10|6.40|
|katsura 11|39.63|
|katsura 12|308.92|
|katsura 13| - |
|noon 7|1.82|
|noon 8|17.75|
|noon 9|225.21|
|noon 10| - |
|reimer 6|0.28|
|reimer 7|7.03|
|reimer 8|272.15|
|reimer 9| - |

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
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
