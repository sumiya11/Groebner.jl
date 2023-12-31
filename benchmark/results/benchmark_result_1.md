## Benchmark results

2023-12-31T02:35:41.951

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|-----|---|---|---|---|---|
|cyclic 7|0.14|0.10|0.15|0.50|1.97|
|cyclic 8|1.54|1.39|1.42|11.09|93.93|
|cyclic 9| - |340.96|244.32| - | - |
|cyclic 10| - | - | - | - | - |
|dummy|0.00|0.01|0.01|0.01|0.00|
|eco 11|0.37|0.40|0.56|3.64|48.27|
|eco 12|2.89|3.07|3.01|17.53| - |
|eco 13|16.31|21.23|20.07|88.54| - |
|eco 14|244.72|231.91|231.99| - | - |
|henrion 5|0.00|0.01|0.02|0.03|0.01|
|henrion 6|0.08|0.05|0.08|0.37|0.40|
|henrion 7|3.13|6.19|3.86|33.34|143.36|
|katsura 10|0.92|1.37|1.22|7.16|180.71|
|katsura 11|10.84|16.80|9.33|44.11| - |
|katsura 12|81.85|129.18|81.17|351.57| - |
|katsura 13| - | - | - | - | - |
|noon 7|0.27|0.21|0.24|3.03|0.48|
|noon 8|2.21|2.09|2.06|19.09|4.84|
|noon 9|22.12|22.13|24.80|289.53|61.24|
|noon 10| - |266.60|339.89| - | - |
|reimer 6|0.05|0.06|0.08|0.32|18.55|
|reimer 7|1.25|1.80|1.24|8.30| - |
|reimer 8|34.03|50.59|46.85|318.63| - |
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
