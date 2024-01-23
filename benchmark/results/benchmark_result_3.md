## Benchmark results

2024-01-23T01:56:04.556

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: The rationals

- Workers: 16
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 9|1.22|4.87|1.50|
|chandra 10|2.92|21.35|4.51|
|chandra 11|10.60|121.88|21.99|
|chandra 12|87.26| - |118.94|
|chandra 13| - | - | - |
|cyclic 7|0.88|1.37|1.37|
|cyclic 8|21.57|21.14|26.78|
|cyclic 9| - | - | - |
|dummy|0.00|0.01|0.03|
|eco 10|0.40|0.66|0.62|
|eco 11|2.88|4.64|3.85|
|eco 12|20.95|31.94|28.75|
|eco 13| - |520.63|380.05|
|henrion 6|0.53|1.69|0.63|
|henrion 7| - | - |428.35|
|ipp|86.51|148.38|20.31|
|katsura 9|2.70|7.76|2.77|
|katsura 10|20.97|76.62|18.49|
|katsura 11| - | - |223.14|
|noon 7|0.48|1.03|1.60|
|noon 8|3.95|8.21|12.65|
|noon 9|26.13|67.71|110.47|
|reimer 6|0.42|0.65|0.31|
|reimer 7|7.46|25.56|6.69|
|reimer 8| - | - |339.74|

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
