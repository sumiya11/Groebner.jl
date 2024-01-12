## Benchmark results

2024-01-12T09:04:15.474

Benchmarked backend: groebner

Benchmark suite: HC

- Workers: 4
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|boon|0.01|
|chandra 5|0.01|
|chandra 6|0.04|
|chandra 7|0.22|
|chandra 8|0.44|
|chandra 9|2.52|
|chandra 10|5.07|
|chandra 11|17.69|
|chandra 12|108.19|
|ipp|185.73|
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
