## Benchmark results

2024-08-16T09:26:06.425

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|BIOMD0000000103|0.18|
|BIOMD0000000123|0.77|
|Cholera|42.38|
|Goodwin (w.)|311.96|
|HIV2|2.00|
|NFkB (w.)|356.92|
|alea6|0.14|
|bayes148|35.39|
|cyclic 7|0.10|
|cyclic 8|1.42|
|cyclic 9|155.86|
|dummy|0.00|
|eco 11|0.55|
|eco 12|2.14|
|eco 13|8.98|
|eco 14|106.58|
|gametwo2|10.59|
|henrion 5|0.00|
|henrion 6|0.05|
|henrion 7|2.22|
|henrion 8| - |
|jason210|5.50|
|katsura 10|0.78|
|katsura 11|5.81|
|katsura 12|33.74|
|katsura 13|723.35|
|mayr42|48.07|
|noon 7|0.15|
|noon 8|1.55|
|noon 9|14.02|
|noon 10|230.04|
|reimer 6|0.05|
|reimer 7|1.57|
|reimer 8|18.06|
|yang1|16.29|

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
