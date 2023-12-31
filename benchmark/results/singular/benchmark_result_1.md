## Benchmark results

2023-12-31T01:56:55.610

Benchmarked backend: singular

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.97|
|cyclic 8|93.93|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.00|
|eco 11|48.27|
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.40|
|henrion 7|143.36|
|katsura 10|180.71|
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.48|
|noon 8|4.84|
|noon 9|61.24|
|noon 10| - |
|reimer 6|18.55|
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
