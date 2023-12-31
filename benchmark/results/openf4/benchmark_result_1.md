## Benchmark results

2023-12-31T02:35:33.899

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.50|
|cyclic 8|11.09|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.01|
|eco 11|3.64|
|eco 12|17.53|
|eco 13|88.54|
|eco 14| - |
|henrion 5|0.03|
|henrion 6|0.37|
|henrion 7|33.34|
|katsura 10|7.16|
|katsura 11|44.11|
|katsura 12|351.57|
|katsura 13| - |
|noon 7|3.03|
|noon 8|19.09|
|noon 9|289.53|
|noon 10| - |
|reimer 6|0.32|
|reimer 7|8.30|
|reimer 8|318.63|
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
