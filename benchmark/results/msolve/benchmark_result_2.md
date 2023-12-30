## Benchmark results

2023-12-29T23:46:11.416

Benchmarked backend: msolve
Benchmark suite: The rationals

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|2.60|
|cyclic 8|68.27|
|dummy|0.03|
|eco 10|0.91|
|eco 11|6.26|
|eco 12|40.62|
|henrion 6|5.71|
|katsura 10|43.81|
|katsura 11|487.08|
|katsura 9|5.33|
|noon 8|18.06|
|noon 9|182.39|

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
