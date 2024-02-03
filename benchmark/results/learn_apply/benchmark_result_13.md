## Benchmark results

2024-02-03T06:11:23.163

Benchmarked backend: learn_apply

Benchmark suite: Learn-apply, other

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|Apply 8x, s|
|:----|---|---|---|---|---|
|BIOMD0000000103|0.24|0.67|0.08|0.13|3.39|
|BIOMD0000000123|0.85|1.03|0.18|0.19|3.48|
|Cholera|46.50|148.29|20.60|37.53|73.88|
|HIV2|2.26|2.64|0.24|0.34|3.70|
|NFkB| - | - | - | - | - |
|Pharm_with_weights| - | - | - | - | - |
|SEAIJRC| - | - | - | - | - |
|SEIRP|76.31|409.59|41.65|80.22|146.57|
|SIWR| - | - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|3.21|
|ipp|0.01|0.01|0.00|0.00|3.21|

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
