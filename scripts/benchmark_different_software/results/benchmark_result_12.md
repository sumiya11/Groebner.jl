## Benchmark results

2024-02-04T14:07:22.902

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: 2^30+3, larger

- Workers: 2
- Timeout: 6000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|dummy|0.00|0.01|0.02|
|eco 15|2115.58|1985.55|2465.36|
|reimer 9|2050.32|4317.16|3354.00|

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
