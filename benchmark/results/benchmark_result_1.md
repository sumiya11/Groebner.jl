## Benchmark results

2024-01-26T15:22:47.295

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 11|0.48|0.76|0.88|
|chandra 12|2.73|2.89|3.84|
|chandra 13|14.95|14.09|20.56|
|chandra 14|98.81|82.15|133.08|
|cyclic 7|0.09|0.24|0.15|
|cyclic 8|1.17|1.38|1.36|
|cyclic 9|120.54|195.46|118.97|
|cyclic 10| - | - | - |
|dummy|0.00|0.18|0.01|
|eco 11|0.30|0.53|0.52|
|eco 12|1.98|2.08|2.59|
|eco 13|9.00|10.52|13.42|
|eco 14|93.24|86.25|127.94|
|henrion 5|0.00|0.20|0.02|
|henrion 6|0.03|0.22|0.06|
|henrion 7|1.97|5.31|3.26|
|henrion 8| - | - | - |
|katsura 10|0.75|1.21|1.08|
|katsura 11|5.14|7.51|6.15|
|katsura 12|35.77|49.28|43.52|
|katsura 13| - |634.70| - |
|noon 7|0.17|0.36|0.22|
|noon 8|1.37|1.75|1.88|
|noon 9|13.44|15.72|16.01|
|noon 10| - |225.13| - |
|noon 11| - | - | - |
|reimer 6|0.05|0.23|0.08|
|reimer 7|0.69|1.41|1.21|
|reimer 8|18.28|26.74|30.17|
|reimer 9| - | - | - |

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
