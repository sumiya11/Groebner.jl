## Benchmark results

2024-01-10T05:06:25.318

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.14|
|cyclic 8|1.40|
|cyclic 9|284.01|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.51|
|eco 12|2.84|
|eco 13|25.01|
|eco 14|275.17|
|henrion 5|0.03|
|henrion 6|0.08|
|henrion 7|3.93|
|katsura 10|1.01|
|katsura 11|9.87|
|katsura 12|110.69|
|katsura 13| - |
|noon 7|0.22|
|noon 8|1.70|
|noon 9|24.81|
|noon 10|362.62|
|reimer 6|0.08|
|reimer 7|1.23|
|reimer 8|49.88|
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
