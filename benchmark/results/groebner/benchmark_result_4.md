## Benchmark results

2024-01-23T15:35:24.929

Benchmarked backend: groebner

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 300 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.00|
|Biohydrogenation_with_weights|0.01|
|COVID_m1|0.64|
|COVID_m1_with_weights|0.17|
|ChemicalReactionNetwork|0.32|
|ChemicalReactionNetwork_with_weights|0.51|
|Cholera|65.06|
|Cholera_with_weights|82.05|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.00|
|DAISY_mamil3|0.00|
|DAISY_mamil3_with_weights|0.00|
|DAISY_mamil4|0.29|
|DAISY_mamil4_with_weights|0.30|
|Goodwin| - |
|Goodwin_with_weights| - |
|HIV|0.03|
|HIV2|3.53|
|HIV2_with_weights|2.51|
|HIV_with_weights|0.01|
|LV|0.00|
|LV_with_weights|0.00|
|Lipolysis|0.00|
|Lipolysis_with_weights|0.00|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose|0.01|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.00|
|SEIR2|0.02|
|SEIR2_with_weights|0.03|
|SEIRP| - |
|SEIRP_with_weights|32.86|
|SEIR_with_weights|0.00|
|SIRSForced|0.87|
|SIRSForced_with_weights|0.25|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.00|
|SlowFast|0.00|
|SlowFast_with_weights|0.00|
|Treatment|0.01|
|Treatment_with_weights|0.01|
|Tumor|0.12|
|Tumor_with_weights|0.17|

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
