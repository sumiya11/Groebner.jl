## Benchmark results

2024-01-10T05:34:04.318

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|-----|---|---|---|---|---|
|cyclic 7|0.12|0.10|0.14|0.53|1.84|
|cyclic 8|1.75|1.35|1.40|13.38|68.34|
|cyclic 9| - |369.73|284.01| - | - |
|cyclic 10| - | - | - | - | - |
|dummy|0.00|0.01|0.01|0.02|0.00|
|eco 11|0.53|0.34|0.51|3.60|37.99|
|eco 12|2.98|2.20|2.84|21.04| - |
|eco 13|23.12|27.60|25.01|136.21| - |
|eco 14| - |262.95|275.17| - | - |
|henrion 5|0.00|0.01|0.03|0.03|0.01|
|henrion 6|0.05|0.06|0.08|0.69|0.39|
|henrion 7|3.78|6.73|3.93|54.26|126.98|
|katsura 10|1.11|1.23|1.01|15.26|120.02|
|katsura 11|12.99|18.69|9.87|59.19| - |
|katsura 12|135.10|165.60|110.69|469.32| - |
|katsura 13| - | - | - | - | - |
|noon 7|0.29|0.19|0.22|3.34|0.45|
|noon 8|2.36|1.69|1.70|41.71|4.65|
|noon 9|25.90|23.04|24.81|482.99|57.79|
|noon 10| - |264.66|362.62| - | - |
|reimer 6|0.17|0.06|0.08|0.47|17.06|
|reimer 7|1.45|1.54|1.23|13.92| - |
|reimer 8|48.70|67.63|49.88|501.20| - |
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
