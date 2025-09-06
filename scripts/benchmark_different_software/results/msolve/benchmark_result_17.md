## Benchmark results

2024-02-07T07:14:53.985

Benchmarked backend: msolve

Benchmark suite: SIAN, QQ

- Workers: 8
- Timeout: 7200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Cholera| - |
|Goodwin (w.)| - |
|HIV2|146.70|
|NFkB (w.)| - |
|Pharm (w.)| - |
|SEIRP| - |
|dummy|0.03|

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
