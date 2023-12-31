## Benchmark results

2023-12-31T14:11:04.046

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7|1.75|
|cyclic 8|30.09|
|dummy|0.02|
|eco 10|0.78|
|eco 11|5.63|
|eco 12|52.84|
|henrion 6| - |
|katsura 9|9.53|
|katsura 10|136.96|
|katsura 11| - |
|noon 8|11.05|
|noon 9|122.91|

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
