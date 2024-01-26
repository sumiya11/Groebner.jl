## Benchmark results

2024-01-22T20:55:01.034

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4"]

Benchmark suite: dummy benchmark set

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|
|:----|---|---|---|---|
|dummy 1|0.00|0.02|0.03| - |
|dummy 2|0.00|0.05|0.03| - |

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
