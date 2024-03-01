## Benchmark results

2024-02-04T04:00:53.343

Benchmarked backend: singular

Benchmark suite: Learn-apply, other

- Workers: 4
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|BIOMD0000000103|19.58|
|BIOMD0000000123|4.12|
|Cholera| - |
|HIV2|201.92|
|NFkB| - |
|Pharm_with_weights| - |
|SEAIJRC| - |
|SEIRP| - |
|SIWR| - |
|dummy|0.00|
|ipp|0.08|

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
