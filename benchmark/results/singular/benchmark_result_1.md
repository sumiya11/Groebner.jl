## Benchmark results

2023-12-30T22:24:29.558

Benchmarked backend: singular

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.98|
|cyclic 8|82.18|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.00|
|eco 11|46.46|
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.40|
|henrion 7|142.64|
|katsura 10|180.66|
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.49|
|noon 8|4.93|
|noon 9|60.44|
|noon 10| - |
|reimer 6|19.50|
|reimer 7| - |
|reimer 8| - |
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
