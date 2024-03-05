## Benchmark results

2024-02-03T04:20:37.617

Benchmarked backend: msolve

Benchmark suite: SIAN, the rationals

- Workers: 10
- Timeout: 3000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.42|
|Biohydrogenation_with_weights|0.36|
|COVID_m1|10.85|
|COVID_m1_with_weights|8.36|
|ChemicalReactionNetwork|5.17|
|ChemicalReactionNetwork_with_weights|3.84|
|Cholera|1841.30|
|Cholera_with_weights|1639.68|
|DAISY_ex3|0.06|
|DAISY_ex3_with_weights|0.10|
|DAISY_mamil3|0.08|
|DAISY_mamil3_with_weights|0.11|
|DAISY_mamil4|2.72|
|DAISY_mamil4_with_weights|7.22|
|Goodwin| - |
|Goodwin_with_weights|2240.06|
|HIV|0.30|
|HIV2|67.04|
|HIV2_with_weights|41.33|
|HIV_with_weights|0.30|
|LV|0.04|
|LV_with_weights|0.04|
|Lipolysis|0.05|
|Lipolysis_with_weights|0.04|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose|0.04|
|OralGlucose_with_weights|0.07|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.06|
|SEIR2|0.16|
|SEIR2_with_weights|0.44|
|SEIRP|2450.64|
|SEIRP_with_weights|313.82|
|SEIR_with_weights|0.07|
|SIRSForced|11.94|
|SIRSForced_with_weights|8.93|
|SIR_R0|0.04|
|SIR_R0_with_weights|0.05|
|SlowFast|0.04|
|SlowFast_with_weights|0.04|
|Treatment|0.14|
|Treatment_with_weights|0.18|
|Tumor|0.82|
|Tumor_with_weights|10.30|

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
