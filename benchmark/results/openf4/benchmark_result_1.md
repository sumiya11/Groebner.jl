## Benchmark results

2024-02-07T07:09:41.171

Benchmarked backend: openf4

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|cyclic 7|0.51|
|cyclic 8|10.14|
|cyclic 9| - |
|dummy|0.02|
|eco 11|3.91|
|eco 12|26.23|
|eco 13|113.64|
|eco 14|1106.28|
|henrion 5|0.04|
|henrion 6|0.52|
|henrion 7|51.99|
|henrion 8| - |
|katsura 10|10.33|
|katsura 11|73.04|
|katsura 12|460.79|
|katsura 13| - |
|noon 7|3.10|
|noon 8|40.31|
|noon 9|413.75|
|reimer 6|0.40|
|reimer 7|13.42|
|reimer 8|439.16|

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
