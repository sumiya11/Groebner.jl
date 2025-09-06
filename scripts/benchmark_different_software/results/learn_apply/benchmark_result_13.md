## Benchmark results

2024-02-09T15:31:14.835

Benchmarked backend: learn_apply

Benchmark suite: Learn-apply, other

- Workers: 8
- Timeout: 10000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|Apply 8x, s|
|:----|---|---|---|---|---|
|BIOMD0000000103|0.24|0.58|0.11|0.12|4.39|
|BIOMD0000000123|1.17|1.25|0.16|0.19|4.11|
|Cholera|107.16|376.82|51.06|117.03|224.11|
|HIV2|2.84|2.89|0.29|0.31|4.18|
|NFkB| - | - | - | - | - |
|Pharm_with_weights| - | - | - | - | - |
|SEIRP|237.51|1307.65|125.60|308.66|597.90|
|dummy|0.00|0.00|0.00|0.00|3.28|
|ipp|0.01|0.01|0.00|0.00|3.40|

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
