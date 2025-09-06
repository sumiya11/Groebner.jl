## Benchmark results

2024-08-16T10:16:31.371

Benchmarked backend: mgb

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|BIOMD0000000103|0.19|
|BIOMD0000000123|1.71|
|Cholera|62.09|
|Goodwin (w.)|292.68|
|HIV2|3.13|
|NFkB (w.)|413.70|
|alea6|0.11|
|bayes148|31.36|
|cyclic 7|0.10|
|cyclic 8|1.32|
|cyclic 9|193.91|
|dummy|0.00|
|eco 11|0.39|
|eco 12|2.30|
|eco 13|14.70|
|eco 14|126.38|
|gametwo2|19.09|
|henrion 5|0.01|
|henrion 6|0.06|
|henrion 7|4.20|
|henrion 8|1624.98|
|jason210|8.23|
|katsura 10|1.17|
|katsura 11|11.35|
|katsura 12|61.22|
|katsura 13|1056.05|
|mayr42|68.16|
|noon 7|0.18|
|noon 8|1.91|
|noon 9|17.26|
|noon 10|200.18|
|reimer 6|0.05|
|reimer 7|1.19|
|reimer 8|37.34|
|yang1|47.20|

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
