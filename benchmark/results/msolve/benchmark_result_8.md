## Benchmark results

2024-01-12T09:08:59.948

Benchmarked backend: msolve

Benchmark suite: HC modulo 2^30 + 3

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 5|0.02|
|chandra 6|0.02|
|chandra 7|0.03|
|chandra 8|0.04|
|chandra 9|0.08|
|chandra 10|0.20|
|chandra 11|0.95|
|chandra 12|5.16|
|chandra 13|33.15|

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
