## Benchmark results

2024-01-12T12:10:34.027

Benchmarked backend: groebner

Benchmark suite: The rationals

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|1.92|
|chandra 10|4.11|
|chandra 11|14.09|
|chandra 12|90.53|
|chandra 13|548.26|
|cyclic 7|0.94|
|cyclic 8|25.17|
|cyclic 9| - |
|dummy|0.00|
|eco 10|0.50|
|eco 11|3.42|
|eco 12|28.89|
|eco 13|359.91|
|henrion 6|0.66|
|henrion 7|394.70|
|ipp|134.31|
|katsura 9|4.11|
|katsura 10|23.64|
|katsura 11|319.80|
|noon 7|0.97|
|noon 8|5.82|
|noon 9|42.68|
|reimer 6|0.54|
|reimer 7|10.80|
|reimer 8|816.37|

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
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
