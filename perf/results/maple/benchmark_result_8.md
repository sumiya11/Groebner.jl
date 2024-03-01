## Benchmark results

2024-01-12T09:07:56.319

Benchmarked backend: maple

Benchmark suite: HC modulo 2^30 + 3

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 5|0.01|
|chandra 6|0.01|
|chandra 7|0.01|
|chandra 8|0.02|
|chandra 9|0.05|
|chandra 10|0.16|
|chandra 11|0.81|
|chandra 12|4.55|
|chandra 13|26.50|

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
