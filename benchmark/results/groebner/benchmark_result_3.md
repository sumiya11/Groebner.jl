## Benchmark results

2024-01-26T16:00:02.075

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.18|
|chandra 10|2.88|
|chandra 11|10.28|
|chandra 12|82.21|
|chandra 13| - |
|cyclic 7|0.80|
|cyclic 8|20.81|
|cyclic 9| - |
|dummy|0.00|
|eco 10|0.44|
|eco 11|3.04|
|eco 12|20.65|
|eco 13|320.03|
|henrion 6|0.53|
|henrion 7| - |
|ipp|61.74|
|katsura 9|2.79|
|katsura 10|20.30|
|katsura 11|297.12|
|noon 7|0.46|
|noon 8|3.84|
|noon 9|24.66|
|reimer 6|0.41|
|reimer 7|6.79|
|reimer 8| - |

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
