## Benchmark results

2024-01-12T13:26:27.874

Benchmarked backend: msolve

Benchmark suite: HC

- Workers: 8
- Timeout: 200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|boon|0.04|
|chandra 5|0.03|
|chandra 6|0.04|
|chandra 7|0.04|
|chandra 8|0.04|
|chandra 9|0.04|
|chandra 10|0.04|
|chandra 11|0.05|
|chandra 12|0.04|
|ipp|0.07|
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
