## Benchmark results

2024-01-10T03:46:27.749

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.12|
|cyclic 8|1.75|
|cyclic 9| - |
|cyclic 10| - |
|dummy|0.00|
|eco 11|0.53|
|eco 12|2.98|
|eco 13|23.12|
|eco 14| - |
|henrion 5|0.00|
|henrion 6|0.05|
|henrion 7|3.78|
|katsura 10|1.11|
|katsura 11|12.99|
|katsura 12|135.10|
|katsura 13| - |
|noon 7|0.29|
|noon 8|2.36|
|noon 9|25.90|
|noon 10| - |
|reimer 6|0.17|
|reimer 7|1.45|
|reimer 8|48.70|
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
