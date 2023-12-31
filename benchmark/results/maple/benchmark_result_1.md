## Benchmark results

2023-12-31T13:37:04.842

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.10|
|cyclic 8|1.38|
|cyclic 9|339.41|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.39|
|eco 12|3.08|
|eco 13|21.31|
|eco 14|233.16|
|henrion 5|0.02|
|henrion 6|0.05|
|henrion 7|6.21|
|katsura 10|1.42|
|katsura 11|16.90|
|katsura 12|129.49|
|katsura 13|1170.66|
|noon 7|0.20|
|noon 8|2.18|
|noon 9|21.70|
|noon 10|268.77|
|reimer 6|0.07|
|reimer 7|1.67|
|reimer 8|50.93|
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
