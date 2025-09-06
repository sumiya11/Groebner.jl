## Benchmark results

2024-01-12T09:24:27.518

Benchmarked backend: singular

Benchmark suite: HC

- Workers: 4
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|boon|0.00|
|chandra 5|0.00|
|chandra 6|0.01|
|chandra 7|0.13|
|chandra 8|0.59|
|chandra 9|3.14|
|chandra 10| - |
|chandra 11| - |
|chandra 12| - |
|ipp| - |
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
