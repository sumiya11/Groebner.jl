## Benchmark results

2024-01-12T13:17:56.773

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: The rationals

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 9|1.92|5.85|WA|
|chandra 10|4.11|27.81|WA|
|chandra 11|14.09|151.60|WA|
|chandra 12|90.53|1166.48|WA|
|chandra 13|548.26| - |WA|
|cyclic 7|0.94|1.50|1.07|
|cyclic 8|25.17|24.92|26.84|
|cyclic 9| - | - | - |
|dummy|0.00|0.01|0.03|
|eco 10|0.50|0.72|0.57|
|eco 11|3.42|4.72|3.67|
|eco 12|28.89|39.68|25.40|
|eco 13|359.91|532.73|443.58|
|henrion 6|0.66|1.96|3.80|
|henrion 7|394.70| - | - |
|ipp|134.31|196.19|WA|
|katsura 9|4.11|9.07|2.77|
|katsura 10|23.64|102.40|20.62|
|katsura 11|319.80|1479.51|244.90|
|noon 7|0.97|1.02|1.56|
|noon 8|5.82|9.62|12.76|
|noon 9|42.68|102.74|102.05|
|reimer 6|0.54|0.65|0.28|
|reimer 7|10.80|37.55|7.06|
|reimer 8|816.37| - |364.48|

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
