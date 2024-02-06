## Benchmark results

2024-02-04T04:59:16.835

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 11|0.59|0.73|0.85|
|chandra 12|3.46|3.35|3.91|
|chandra 13|18.84|17.89|20.44|
|cyclic 7|0.11|0.11|0.14|
|cyclic 8|1.30|1.41|1.38|
|cyclic 9|150.87|270.81|108.98|
|dummy|0.00|0.01|0.01|
|eco 11|0.35|0.36|0.57|
|eco 12|2.05|2.23|2.57|
|eco 13|11.75|14.31|12.64|
|eco 14|133.04|140.89|129.87|
|henrion 5|0.00|0.01|0.02|
|henrion 6|0.03|0.05|0.08|
|henrion 7|2.34|3.92|3.04|
|henrion 8|871.57|1190.94|563.32|
|katsura 10|0.82|1.19|1.07|
|katsura 11|6.54|10.79|6.22|
|katsura 12|50.80|69.45|39.49|
|katsura 13|483.72|922.76|267.54|
|katsura 14| - | - | - |
|noon 7|0.19|0.19|0.25|
|noon 8|1.72|1.97|1.73|
|noon 9|16.14|19.40|16.28|
|noon 10|188.68|209.29|170.61|
|reimer 6|0.06|0.06|0.07|
|reimer 7|1.02|1.47|1.12|
|reimer 8|25.37|47.17|28.36|

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
