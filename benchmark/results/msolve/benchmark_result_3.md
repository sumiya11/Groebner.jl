## Benchmark results

2023-12-31T14:45:53.543

Benchmarked backend: msolve

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.20|
|cyclic 8|25.64|
|dummy|0.03|
|eco 10|0.62|
|eco 11|3.58|
|eco 12|27.44|
|henrion 6|3.95|
|katsura 9|2.35|
|katsura 10|17.89|
|katsura 11|197.02|
|noon 8|11.92|
|noon 9|126.53|

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
