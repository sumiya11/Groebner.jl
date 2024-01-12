## Benchmark results

2024-01-12T09:17:58.879

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: HC modulo 2^30 + 3

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|:----|---|---|---|---|---|
|chandra 5|0.00|0.01|0.02|0.01|0.00|
|chandra 6|0.00|0.01|0.02|0.03|0.00|
|chandra 7|0.00|0.01|0.03|0.03|0.01|
|chandra 8|0.01|0.02|0.04|0.08|0.06|
|chandra 9|0.04|0.05|0.08|0.56|0.40|
|chandra 10|0.32|0.16|0.20|1.70|2.62|
|chandra 11|0.87|0.81|0.95|9.92|23.09|
|chandra 12|4.65|4.55|5.16|71.29|183.02|
|chandra 13|32.91|26.50|33.15|503.61| - |

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
