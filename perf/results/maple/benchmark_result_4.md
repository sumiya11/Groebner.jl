## Benchmark results

2024-02-02T22:44:27.037

Benchmarked backend: maple

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 10
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.18|
|Biohydrogenation_with_weights|0.25|
|COVID_m1|0.82|
|COVID_m1_with_weights|0.28|
|ChemicalReactionNetwork|0.41|
|ChemicalReactionNetwork_with_weights|0.37|
|Cholera|40.30|
|Cholera_with_weights|44.84|
|DAISY_ex3|0.14|
|DAISY_ex3_with_weights|0.16|
|DAISY_mamil3|0.17|
|DAISY_mamil3_with_weights|0.16|
|DAISY_mamil4|0.56|
|DAISY_mamil4_with_weights|0.57|
|Goodwin| - |
|Goodwin_with_weights|206.42|
|HIV|0.18|
|HIV2|3.18|
|HIV2_with_weights|1.57|
|HIV_with_weights|0.17|
|LV|0.16|
|LV_with_weights|0.18|
|Lipolysis|0.15|
|Lipolysis_with_weights|0.20|
|NFkB|878.10|
|NFkB_with_weights|249.35|
|OralGlucose|0.17|
|OralGlucose_with_weights|0.18|
|Pharm| - |
|Pharm_with_weights|1310.25|
|SEIR|0.15|
|SEIR2|0.22|
|SEIR2_with_weights|0.35|
|SEIRP|103.72|
|SEIRP_with_weights|38.65|
|SEIR_with_weights|0.20|
|SIRSForced|1.11|
|SIRSForced_with_weights|0.31|
|SIR_R0|0.23|
|SIR_R0_with_weights|0.21|
|SlowFast|0.19|
|SlowFast_with_weights|0.21|
|Treatment|0.22|
|Treatment_with_weights|0.19|
|Tumor|0.25|
|Tumor_with_weights|0.23|

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
