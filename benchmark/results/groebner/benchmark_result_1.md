## Benchmark results

2024-01-26T15:00:05.266

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.48|
|chandra 12|2.73|
|chandra 13|14.95|
|chandra 14|98.81|
|cyclic 7|0.09|
|cyclic 8|1.17|
|cyclic 9|120.54|
|cyclic 10| - |
|dummy|0.00|
|eco 11|0.30|
|eco 12|1.98|
|eco 13|9.00|
|eco 14|93.24|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|1.97|
|henrion 8| - |
|katsura 10|0.75|
|katsura 11|5.14|
|katsura 12|35.77|
|katsura 13| - |
|noon 7|0.17|
|noon 8|1.37|
|noon 9|13.44|
|noon 10| - |
|noon 11| - |
|reimer 6|0.05|
|reimer 7|0.69|
|reimer 8|18.28|
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
