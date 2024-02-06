## Benchmark results

2024-02-03T03:33:25.755

Benchmarked backend: learn_apply

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|F4, s|Learn, s|Apply, s|Apply 4x, s|Apply 8x, s|
|:----|---|---|---|---|---|
|chandra 11|0.51|1.19|0.08|0.20|3.59|
|chandra 12|2.91|7.91|0.26|0.61|4.49|
|chandra 13|15.21|50.77|0.87|2.41|8.11|
|chandra 14| - | - | - | - | - |
|cyclic 7|0.11|0.21|0.03|0.06|3.09|
|cyclic 8|1.35|4.04|0.60|0.96|4.71|
|cyclic 9| - | - | - | - | - |
|dummy|0.00|0.00|0.00|0.00|3.09|
|eco 11|0.36|1.19|0.16|0.23|3.53|
|eco 12|2.13|11.10|0.98|1.59|6.23|
|eco 13|9.18|85.75|3.39|5.52|13.53|
|eco 14| - | - | - | - | - |
|henrion 5|0.00|0.00|0.00|0.00|3.07|
|henrion 6|0.03|0.08|0.01|0.03|3.11|
|henrion 7|2.01|5.55|0.82|1.45|5.65|
|henrion 8| - | - | - | - | - |
|katsura 10|0.86|5.19|0.33|0.63|4.33|
|katsura 11|5.16|50.12|2.27|3.98|10.99|
|katsura 12| - | - | - | - | - |
|katsura 13| - | - | - | - | - |
|noon 7|0.16|0.21|0.05|0.05|3.32|
|noon 8|1.32|1.53|0.20|0.39|3.94|
|noon 9|13.21|13.97|1.40|2.99|8.90|
|noon 10| - | - | - | - | - |
|reimer 6|0.07|0.06|0.01|0.02|3.10|
|reimer 7|0.78|1.34|0.12|0.18|3.49|
|reimer 8|18.82|44.66|1.96|2.97|8.77|

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
