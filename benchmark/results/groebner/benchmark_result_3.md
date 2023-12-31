## Benchmark results

2023-12-31T14:39:14.826

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.19|
|cyclic 8|27.54|
|dummy|0.00|
|eco 10|0.68|
|eco 11|3.53|
|eco 12|31.14|
|henrion 6|0.76|
|katsura 9|3.73|
|katsura 10|27.37|
|katsura 11|355.52|
|noon 8|5.29|
|noon 9|56.26|

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
