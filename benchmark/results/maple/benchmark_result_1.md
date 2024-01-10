## Benchmark results

2024-01-10T04:46:42.469

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.10|
|cyclic 8|1.35|
|cyclic 9|369.73|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.34|
|eco 12|2.20|
|eco 13|27.60|
|eco 14|262.95|
|henrion 5|0.01|
|henrion 6|0.06|
|henrion 7|6.73|
|katsura 10|1.23|
|katsura 11|18.69|
|katsura 12|165.60|
|katsura 13| - |
|noon 7|0.19|
|noon 8|1.69|
|noon 9|23.04|
|noon 10|264.66|
|reimer 6|0.06|
|reimer 7|1.54|
|reimer 8|67.63|
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
