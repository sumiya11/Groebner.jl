## Benchmark results

2023-12-30T01:46:20.918

Benchmarked backend: groebner
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.13|
|cyclic 8|1.38|
|cyclic 9|219.54|
|dummy|0.00|
|eco 11|0.44|
|eco 12|2.78|
|eco 13|17.33|
|eco 14|195.25|
|henrion 5|0.00|
|henrion 6|0.07|
|henrion 7|3.01|
|katsura 10|0.96|
|katsura 11|8.85|
|katsura 12|72.87|
|katsura 13|601.87|
|noon 10|291.45|
|noon 7|0.23|
|noon 8|1.96|
|noon 9|23.12|
|reimer 6|0.05|
|reimer 7|1.34|
|reimer 8|37.38|
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
