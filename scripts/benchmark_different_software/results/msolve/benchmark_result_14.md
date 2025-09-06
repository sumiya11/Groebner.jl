## Benchmark results

2024-02-03T14:14:14.566

Benchmarked backend: msolve

Benchmark suite: Other, modulo 2^30 + 3

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|alea6|0.16|
|bayes148|60.78|
|dummy|0.01|
|gametwo2|20.23|
|jason210|9.17|
|mayr42|90.57|
|yang1|69.93|

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
