## Benchmark results

2024-02-07T06:42:34.534

Benchmarked backend: maple

Benchmark suite: SIAN, QQ

- Workers: 8
- Timeout: 7200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Cholera|527.64|
|Goodwin (w.)|2120.47|
|HIV2|23.85|
|NFkB (w.)|3411.72|
|Pharm (w.)| - |
|SEIRP|1235.47|
|dummy|0.01|

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
