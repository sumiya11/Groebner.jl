## Benchmark results

2024-02-02T16:46:36.124

Benchmarked backend: groebner

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 10
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.00|
|Biohydrogenation_with_weights|0.00|
|COVID_m1|0.50|
|COVID_m1_with_weights|0.12|
|ChemicalReactionNetwork|0.23|
|ChemicalReactionNetwork_with_weights|0.24|
|Cholera|43.00|
|Cholera_with_weights|45.60|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.00|
|DAISY_mamil3|0.02|
|DAISY_mamil3_with_weights|0.00|
|DAISY_mamil4|0.22|
|DAISY_mamil4_with_weights|0.23|
|Goodwin| - |
|Goodwin_with_weights|152.54|
|HIV|0.01|
|HIV2|2.27|
|HIV2_with_weights|1.22|
|HIV_with_weights|0.00|
|LV|0.00|
|LV_with_weights|0.02|
|Lipolysis|0.00|
|Lipolysis_with_weights|0.00|
|NFkB|690.76|
|NFkB_with_weights|171.85|
|OralGlucose|0.00|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights|742.05|
|SEIR|0.03|
|SEIR2|0.00|
|SEIR2_with_weights|0.02|
|SEIRP|76.15|
|SEIRP_with_weights|25.77|
|SEIR_with_weights|0.00|
|SIRSForced|0.63|
|SIRSForced_with_weights|0.13|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.00|
|SlowFast|0.00|
|SlowFast_with_weights|0.00|
|Treatment|0.00|
|Treatment_with_weights|0.01|
|Tumor|0.04|
|Tumor_with_weights|0.06|

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
