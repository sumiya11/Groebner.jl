## Benchmark results

2023-12-30T22:35:22.682

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.10|
|cyclic 8|1.32|
|cyclic 9|333.67|
|cyclic 10| - |
|dummy|0.02|
|eco 11|0.41|
|eco 12|2.46|
|eco 13|21.11|
|eco 14|228.48|
|henrion 5|0.02|
|henrion 6|0.06|
|henrion 7|6.39|
|katsura 10|1.32|
|katsura 11|16.77|
|katsura 12|131.91|
|katsura 13| - |
|noon 7|0.20|
|noon 8|1.90|
|noon 9|21.58|
|noon 10|268.01|
|reimer 6|0.06|
|reimer 7|1.56|
|reimer 8|49.79|
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
