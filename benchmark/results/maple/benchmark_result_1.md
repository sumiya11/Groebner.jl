## Benchmark results

2023-12-30T02:55:08.019

Benchmarked backend: maple
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.10|
|cyclic 8|1.44|
|cyclic 9|297.84|
|dummy|0.01|
|eco 11|0.37|
|eco 12|2.58|
|eco 13|21.81|
|eco 14|203.63|
|henrion 5|0.02|
|henrion 6|0.05|
|henrion 7|5.79|
|katsura 10|1.37|
|katsura 11|17.56|
|katsura 12|135.77|
|katsura 13|1012.59|
|noon 10|255.76|
|noon 7|0.20|
|noon 8|1.95|
|noon 9|20.17|
|reimer 6|0.06|
|reimer 7|1.53|
|reimer 8|57.62|
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
