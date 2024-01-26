## Benchmark results

2024-01-26T16:58:14.560

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: The rationals

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 9|1.18|4.84|1.25|
|chandra 10|2.88|21.17|4.53|
|chandra 11|10.28|111.70|19.83|
|chandra 12|82.21|886.62|105.81|
|chandra 13| - | - | - |
|cyclic 7|0.80|1.43|1.20|
|cyclic 8|20.81|20.50|27.28|
|cyclic 9| - | - | - |
|dummy|0.00|0.14|0.02|
|eco 10|0.44|0.93|0.65|
|eco 11|3.04|4.79|4.02|
|eco 12|20.65|31.59|28.25|
|eco 13|320.03|352.97|314.83|
|henrion 6|0.53|2.03|0.65|
|henrion 7| - | - |404.00|
|ipp|61.74|91.40|14.18|
|katsura 9|2.79|7.70|2.21|
|katsura 10|20.30|76.88|17.29|
|katsura 11|297.12|1218.94|206.39|
|noon 7|0.46|1.09|1.53|
|noon 8|3.84|8.15|11.80|
|noon 9|24.66|66.52|100.94|
|reimer 6|0.41|0.77|0.30|
|reimer 7|6.79|24.22|6.14|
|reimer 8| - | - |358.83|

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
