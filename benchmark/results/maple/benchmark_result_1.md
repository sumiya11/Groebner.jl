## Benchmark results

2024-01-23T15:00:09.314

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.62|
|chandra 12|2.72|
|chandra 13|14.51|
|chandra 14|103.12|
|cyclic 7|0.10|
|cyclic 8|1.20|
|cyclic 9|209.38|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.36|
|eco 12|1.95|
|eco 13|10.94|
|eco 14|125.95|
|henrion 5|0.01|
|henrion 6|0.05|
|henrion 7|3.03|
|henrion 8| - |
|katsura 10|1.26|
|katsura 11|7.33|
|katsura 12|57.21|
|katsura 13|987.88|
|noon 7|0.19|
|noon 8|2.12|
|noon 9|15.12|
|noon 10|238.63|
|noon 11| - |
|reimer 6|0.06|
|reimer 7|1.67|
|reimer 8|28.65|
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
