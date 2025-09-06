## Benchmark results

2024-02-02T17:56:44.407

Benchmarked backend: msolve

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 10
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.02|
|Biohydrogenation_with_weights|0.01|
|COVID_m1|1.72|
|COVID_m1_with_weights|0.30|
|ChemicalReactionNetwork|0.55|
|ChemicalReactionNetwork_with_weights| - |
|Cholera|111.14|
|Cholera_with_weights|102.71|
|DAISY_ex3|0.02|
|DAISY_ex3_with_weights|12.72|
|DAISY_mamil3|0.02|
|DAISY_mamil3_with_weights|1.83|
|DAISY_mamil4|0.61|
|DAISY_mamil4_with_weights| - |
|Goodwin| - |
|Goodwin_with_weights|458.67|
|HIV|0.02|
|HIV2|11.74|
|HIV2_with_weights|3.27|
|HIV_with_weights|0.02|
|LV|0.01|
|LV_with_weights|0.01|
|Lipolysis|0.02|
|Lipolysis_with_weights|0.02|
|NFkB|1414.47|
|NFkB_with_weights| - |
|OralGlucose|0.02|
|OralGlucose_with_weights|0.02|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.02|
|SEIR2|0.03|
|SEIR2_with_weights| - |
|SEIRP|160.54|
|SEIRP_with_weights|43.71|
|SEIR_with_weights|0.07|
|SIRSForced|2.36|
|SIRSForced_with_weights|0.35|
|SIR_R0|0.02|
|SIR_R0_with_weights|0.02|
|SlowFast|0.02|
|SlowFast_with_weights|0.02|
|Treatment|0.03|
|Treatment_with_weights|0.03|
|Tumor|0.11|
|Tumor_with_weights|0.08|

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
