## Benchmark results

2024-01-12T09:17:56.545

Benchmarked backend: openf4

Benchmark suite: HC modulo 2^30 + 3

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 5|0.01|
|chandra 6|0.03|
|chandra 7|0.03|
|chandra 8|0.08|
|chandra 9|0.56|
|chandra 10|1.70|
|chandra 11|9.92|
|chandra 12|71.29|
|chandra 13|503.61|

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
