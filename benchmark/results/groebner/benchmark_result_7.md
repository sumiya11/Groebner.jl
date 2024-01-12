## Benchmark results

2024-01-12T13:17:47.010

Benchmarked backend: groebner

Benchmark suite: HC

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|boon|0.00|
|chandra 5|0.01|
|chandra 6|0.02|
|chandra 7|0.07|
|chandra 8|0.30|
|chandra 9|1.26|
|chandra 10|4.70|
|chandra 11|16.81|
|chandra 12|107.68|
|ipp|195.66|
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
