## Benchmark results

2023-12-30T20:42:13.817

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 20 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|0.11|0.10|0.15|0.51|1.99|
|cyclic 8|1.25|1.19|1.89|9.94| - |
|cyclic 9| - | - | - | - | - |
|cyclic 10| - | - | - | - | - |
|dummy|0.00|0.01|0.01|0.01|0.00|
|eco 11|0.34|0.35|0.80|2.40| - |
|eco 12|2.14|1.82|3.04| - | - |
|eco 13| - |10.02|18.31| - | - |
|eco 14| - | - | - | - | - |
|henrion 5|0.00|0.01|0.02| - |0.01|
|henrion 6|0.05|0.05|0.08| - |0.41|
|henrion 7|2.08|2.76|3.45| - | - |
|katsura 10|0.90|0.98|1.20|10.06| - |
|katsura 11|5.69|6.20|6.51| - | - |
|katsura 12| - | - | - | - | - |
|katsura 13| - | - | - | - | - |
|noon 7|0.17|0.19|0.26|3.25|0.48|
|noon 8|1.94|1.55|1.82| - | - |
|noon 9| - |13.79|16.96| - | - |
|noon 10| - | - | - | - | - |
|reimer 6|0.04|0.06|0.08|0.50| - |
|reimer 7|0.83|0.94|1.12|8.30| - |
|reimer 8| - |24.49| - | - | - |
|reimer 9| - | - | - | - | - |

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
