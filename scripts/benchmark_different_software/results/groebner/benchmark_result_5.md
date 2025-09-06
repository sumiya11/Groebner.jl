## Benchmark results

2024-01-12T03:37:20.131

Benchmarked backend: groebner

Benchmark suite: MQ

- Workers: 4
- Timeout: 20000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|mq_n10_m7_p2_s0|0.43|
|mq_n10_m20_p2_s0|0.01|
|mq_n10_m20_p2_s1|0.01|
|mq_n15_m10_p2_s0|185.03|
|mq_n15_m30_p2_s0|0.63|
|mq_n24_m16_p31_s0| - |
|mq_n34_m68_p31_s0| - |

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
