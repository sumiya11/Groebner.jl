## Benchmark results

2024-08-16T09:48:56.801

Benchmarked backend: maplefgb

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|BIOMD0000000103|0.20|
|BIOMD0000000123|1.77|
|Cholera|53.67|
|Goodwin (w.)|273.11|
|HIV2|3.00|
|NFkB (w.)|434.36|
|alea6|0.11|
|bayes148|30.44|
|cyclic 7|0.10|
|cyclic 8|1.29|
|cyclic 9|195.90|
|dummy|0.02|
|eco 11|0.36|
|eco 12|2.05|
|eco 13|11.58|
|eco 14|104.79|
|gametwo2|13.40|
|henrion 5|0.01|
|henrion 6|0.05|
|henrion 7|3.25|
|henrion 8| - |
|jason210|7.76|
|katsura 10|1.18|
|katsura 11|11.36|
|katsura 12|53.16|
|katsura 13|879.25|
|mayr42|52.03|
|noon 7|0.19|
|noon 8|1.48|
|noon 9|14.80|
|noon 10|182.06|
|reimer 6|0.05|
|reimer 7|1.18|
|reimer 8|31.03|
|yang1|40.33|

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
