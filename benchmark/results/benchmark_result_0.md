## Benchmark results

2024-01-27T13:24:45.699

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: dummy benchmark set

- Workers: 1
- Timeout: 1800 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|:----|---|---|---|---|---|
|dummy 1|0.00|0.22|0.01|0.01|0.00|
|dummy 2|0.00|0.21|0.01|0.01|0.00|

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
* HostCPUFeatures : 0.1.16
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
