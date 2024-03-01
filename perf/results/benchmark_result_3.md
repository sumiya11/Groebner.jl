## Benchmark results

2024-02-06T12:27:11.354

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: The rationals

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|chandra 9|1.03|5.37|1.29|
|chandra 10|3.98|25.70|4.55|
|chandra 11|15.26|129.58|21.31|
|chandra 12|95.65|1052.14|104.74|
|chandra 13|527.64| - |553.96|
|cyclic 7|0.95|1.41|1.13|
|cyclic 8|19.67|23.82|26.10|
|dummy|0.00|0.01|0.02|
|eco 10|0.48|0.71|0.59|
|eco 11|2.69|4.87|3.72|
|eco 12|24.53|35.06|29.19|
|eco 13|288.08|496.48|269.46|
|henrion 6|0.47|1.74|0.67|
|henrion 7|347.42| - |391.25|
|ipp|34.76|173.43|16.87|
|katsura 9|3.17|9.00|2.27|
|katsura 10|18.65|84.81|17.54|
|katsura 11|306.59|1318.31|167.57|
|noon 7|0.81|1.03|1.51|
|noon 8|4.90|9.23|11.37|
|noon 9|30.11|98.78|91.26|
|reimer 6|0.47|0.57|0.32|
|reimer 7|9.08|32.33|6.68|
|reimer 8|838.21| - |256.74|

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
