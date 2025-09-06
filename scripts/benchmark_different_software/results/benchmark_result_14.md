## Benchmark results

2024-02-03T14:35:19.895

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: Other, modulo 2^30 + 3

- Workers: 8
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|:----|---|---|---|---|---|
|alea6|0.22|0.12|0.16|0.62|4.69|
|bayes148|80.74|44.87|60.78| - |7.48|
|dummy|0.00|0.01|0.01|0.01|0.00|
|gametwo2|23.63|17.22|20.23|38.92| - |
|jason210|7.47|8.34|9.17| - |2.39|
|mayr42|107.78|97.52|90.57| - |178.37|
|yang1|27.90|40.42|69.93| - |78.44|

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
