## Benchmark results

2023-12-29T22:11:44.031

Benchmarked backend: maple
Benchmark suite: The rationals

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.74|
|cyclic 8|30.75|
|dummy|0.01|
|eco 10|0.74|
|eco 11|5.58|
|henrion 6| - |
|katsura 10|126.77|
|katsura 9|9.77|
|noon 8|9.73|
|noon 9|122.21|

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
