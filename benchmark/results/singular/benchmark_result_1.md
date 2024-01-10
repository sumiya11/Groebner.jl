## Benchmark results

2024-01-10T04:27:42.467

Benchmarked backend: singular

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.84|
|cyclic 8|68.34|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.00|
|eco 11|37.99|
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5|0.01|
|henrion 6|0.39|
|henrion 7|126.98|
|katsura 10|120.02|
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.45|
|noon 8|4.65|
|noon 9|57.79|
|noon 10| - |
|reimer 6|17.06|
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
