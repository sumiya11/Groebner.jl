## Benchmark results

2024-08-16T10:21:30.684

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|BIOMD0000000103| - |
|BIOMD0000000123| - |
|Cholera| - |
|Goodwin (w.)| - |
|HIV2| - |
|NFkB (w.)| - |
|alea6| - |
|bayes148| - |
|cyclic 7| - |
|cyclic 8| - |
|cyclic 9| - |
|dummy| - |
|eco 11| - |
|eco 12| - |
|eco 13| - |
|eco 14| - |
|gametwo2| - |
|henrion 5| - |
|henrion 6| - |
|henrion 7| - |
|henrion 8| - |
|jason210| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|mayr42| - |
|noon 7| - |
|noon 8| - |
|noon 9| - |
|noon 10| - |
|reimer 6| - |
|reimer 7| - |
|reimer 8| - |
|yang1| - |

*Benchmarking environment:*

* Total RAM (GiB): 188
* Processor: Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
* Julia version: 1.10.2

Versions of the dependencies:

* Primes : 0.5.6
* TimerOutputs : 0.5.24
* PrecompileTools : 1.2.1
* MultivariatePolynomials : 0.5.6
* Combinatorics : 1.0.2
* HostCPUFeatures : 0.1.17
* AbstractAlgebra : 0.42.0
* Nemo : 0.46.0
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
