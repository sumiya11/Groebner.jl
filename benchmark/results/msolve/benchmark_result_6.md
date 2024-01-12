## Benchmark results

2024-01-11T14:21:00.687

Benchmarked backend: msolve

Benchmark suite: SI

- Workers: 4
- Timeout: 500 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|-----|---|
|SEAIJRC|45.93|
|SIWR|5.08|

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
