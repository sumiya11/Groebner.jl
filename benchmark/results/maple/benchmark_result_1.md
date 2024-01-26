## Benchmark results

2024-01-26T15:11:24.864

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.76|
|chandra 12|2.89|
|chandra 13|14.09|
|chandra 14|82.15|
|cyclic 7|0.24|
|cyclic 8|1.38|
|cyclic 9|195.46|
|cyclic 10| - |
|dummy|0.18|
|eco 11|0.53|
|eco 12|2.08|
|eco 13|10.52|
|eco 14|86.25|
|henrion 5|0.20|
|henrion 6|0.22|
|henrion 7|5.31|
|henrion 8| - |
|katsura 10|1.21|
|katsura 11|7.51|
|katsura 12|49.28|
|katsura 13|634.70|
|noon 7|0.36|
|noon 8|1.75|
|noon 9|15.72|
|noon 10|225.13|
|noon 11| - |
|reimer 6|0.23|
|reimer 7|1.41|
|reimer 8|26.74|
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
