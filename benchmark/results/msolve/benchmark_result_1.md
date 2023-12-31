## Benchmark results

2023-12-31T16:06:03.144

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.19|
|cyclic 8|2.52|
|cyclic 9|806.99|
|cyclic 10| - |
|dummy|0.01|
|eco 11|1.12|
|eco 12|8.95|
|eco 13|109.62|
|eco 14|1606.64|
|henrion 5|0.02|
|henrion 6|0.09|
|henrion 7|6.98|
|katsura 10|3.82|
|katsura 11|68.06|
|katsura 12|652.83|
|katsura 13| - |
|noon 7|0.27|
|noon 8|2.20|
|noon 9|25.14|
|noon 10|234.47|
|reimer 6|0.08|
|reimer 7|1.63|
|reimer 8|85.49|
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
