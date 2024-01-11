## Benchmark results

2024-01-11T10:13:26.466

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|
|-----|---|---|---|---|
|cyclic 7|0.08|0.21|0.03|0.05|
|cyclic 8|1.26|4.13|0.55|0.77|
|cyclic 9| - | - | - | - |
|cyclic 10| - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|
|eco 11|0.37|1.16|0.19|0.27|
|eco 12|2.80|12.12|1.11|1.89|
|eco 13|14.62|106.86|4.41|6.30|
|eco 14| - | - | - | - |
|henrion 5|0.00|0.00|0.00|0.00|
|henrion 6|0.08|0.06|0.02|0.06|
|henrion 7|2.67|7.92|1.11|1.70|
|henrion 8| - | - | - | - |
|katsura 10|0.92|5.34|0.35|0.74|
|katsura 11|9.96|113.54|3.65|6.17|
|katsura 12| - | - | - | - |
|katsura 13| - | - | - | - |
|noon 7|0.17|0.30|0.04|0.08|
|noon 8|1.99|2.02|0.33|0.58|
|noon 9|20.80|19.38|1.74|3.55|
|noon 10|164.54|139.06|10.37|25.71|
|reimer 6|0.06|0.08|0.07|0.02|
|reimer 7|1.19|1.74|0.11|0.28|
|reimer 8|30.14|84.51|2.60|4.25|
|reimer 9| - | - | - | - |

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
