## Benchmark results

2023-12-31T02:52:27.789

Benchmarked backends: Any["groebner", "maple", "msolve", "singular"]

Benchmark suite: The rationals

- Workers: 8
- Timeout: 100 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|singular|
|-----|---|---|---|---|
|cyclic 7|1.23|1.69|2.56| - |
|cyclic 8|28.93|33.62|67.23| - |
|dummy|0.00|0.01|0.02|0.00|
|eco 10|0.67|0.92|0.80|30.60|
|eco 11|3.55|5.64|6.26| - |
|eco 12|29.99|52.92|41.68| - |
|henrion 6|0.72| - |6.34|5.22|
|katsura 9|4.06|10.30|5.74| - |
|katsura 10|29.28| - |43.10| - |
|katsura 11| - | - | - | - |
|noon 8|5.31|10.19|17.45|13.60|
|noon 9| - | - | - | - |

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
