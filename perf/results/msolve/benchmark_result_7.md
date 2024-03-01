## Benchmark results

2024-01-13T09:50:12.191

Benchmarked backend: msolve

Benchmark suite: HC

- Workers: 2
- Timeout: 120 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|boon|0.04|
|chandra 2|0.04|
|chandra 3|0.04|
|chandra 4|0.04|
|chandra 5|0.05|
|chandra 6|0.07|
|chandra 7|0.14|
|chandra 8|0.34|
|chandra 9|1.49|
|chandra 10|6.85|
|chandra 11|31.64|
|chandra 12| - |
|ipp|14.48|
|rps10| - |

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
