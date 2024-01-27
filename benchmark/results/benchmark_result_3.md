## Benchmark results

2024-01-27T08:58:45.352

Benchmarked backends: Any["groebner", "maple", "msolve", "singular"]

Benchmark suite: The rationals

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|singular|
|:----|---|---|---|---|
|chandra 9|1.28|5.53|1.44|3.10|
|chandra 10|4.20|25.24|5.25|35.13|
|chandra 11|13.27|137.17|21.87|284.97|
|chandra 12|82.00|927.56|115.73| - |
|chandra 13|415.78| - |709.91| - |
|cyclic 7|1.11|1.60|1.11| - |
|cyclic 8|27.19|22.57|27.03| - |
|cyclic 9| - | - | - | - |
|dummy|0.00|0.16|0.03|0.00|
|eco 10|0.47|0.78|0.66|43.90|
|eco 11|3.16|5.30|4.00| - |
|eco 12|26.92|39.84|29.09| - |
|eco 13|281.37|420.79|344.04| - |
|henrion 6|0.63|1.91|0.63|6.16|
|henrion 7|346.59|3205.31|364.54| - |
|ipp|93.93|82.50|16.44|870.29|
|katsura 9|3.49|8.71|2.52|284.81|
|katsura 10|24.64|93.27|19.68| - |
|katsura 11|269.47|1439.72|198.94| - |
|noon 7|0.66|1.18|1.75|1.11|
|noon 8|5.37|9.32|11.94|13.23|
|noon 9|36.32|75.82|102.60|172.35|
|reimer 6|0.55|0.80|0.28| - |
|reimer 7|10.48|31.54|6.43| - |
|reimer 8|661.18|3425.27|296.64| - |

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
* HostCPUFeatures : 0.1.16
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
