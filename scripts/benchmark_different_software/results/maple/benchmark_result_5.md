## Benchmark results

2024-01-11T13:23:04.731

Benchmarked backend: maple

Benchmark suite: MQ

- Workers: 4
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|-----|---|
|mq_n10_m7_p2_s0|0.35|
|mq_n10_m20_p2_s0|0.02|
|mq_n10_m20_p2_s1|0.03|
|mq_n15_m10_p2_s0| - |
|mq_n15_m30_p2_s0|0.52|
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
