## Benchmark results

2024-08-16T10:21:32.201

Benchmarked backends: Any["groebner", "maplefgb", "mgb", "msolve", "openf4", "singular"]

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maplefgb|mgb|msolve|openf4|singular|
|:----|---|---|---|---|---|---|
|BIOMD0000000103|0.18|0.20|0.19| - | - | - |
|BIOMD0000000123|0.77|1.77|1.71| - | - | - |
|Cholera|42.38|53.67|62.09| - | - | - |
|Goodwin (w.)|311.96|273.11|292.68| - | - | - |
|HIV2|2.00|3.00|3.13| - | - | - |
|NFkB (w.)|356.92|434.36|413.70| - | - | - |
|alea6|0.14|0.11|0.11| - | - | - |
|bayes148|35.39|30.44|31.36| - | - | - |
|cyclic 7|0.10|0.10|0.10| - | - | - |
|cyclic 8|1.42|1.29|1.32| - | - | - |
|cyclic 9|155.86|195.90|193.91| - | - | - |
|dummy|0.00|0.02|0.00| - | - | - |
|eco 11|0.55|0.36|0.39| - | - | - |
|eco 12|2.14|2.05|2.30| - | - | - |
|eco 13|8.98|11.58|14.70| - | - | - |
|eco 14|106.58|104.79|126.38| - | - | - |
|gametwo2|10.59|13.40|19.09| - | - | - |
|henrion 5|0.00|0.01|0.01| - | - | - |
|henrion 6|0.05|0.05|0.06| - | - | - |
|henrion 7|2.22|3.25|4.20| - | - | - |
|henrion 8| - | - |1624.98| - | - | - |
|jason210|5.50|7.76|8.23| - | - | - |
|katsura 10|0.78|1.18|1.17| - | - | - |
|katsura 11|5.81|11.36|11.35| - | - | - |
|katsura 12|33.74|53.16|61.22| - | - | - |
|katsura 13|723.35|879.25|1056.05| - | - | - |
|mayr42|48.07|52.03|68.16| - | - | - |
|noon 7|0.15|0.19|0.18| - | - | - |
|noon 8|1.55|1.48|1.91| - | - | - |
|noon 9|14.02|14.80|17.26| - | - | - |
|noon 10|230.04|182.06|200.18| - | - | - |
|reimer 6|0.05|0.05|0.05| - | - | - |
|reimer 7|1.57|1.18|1.19| - | - | - |
|reimer 8|18.06|31.03|37.34| - | - | - |
|yang1|16.29|40.33|47.20| - | - | - |

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
