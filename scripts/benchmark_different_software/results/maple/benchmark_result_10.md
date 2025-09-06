## Benchmark results

2024-02-03T03:29:11.091

Benchmarked backend: maple

Benchmark suite: SIAN, the rationals

- Workers: 10
- Timeout: 3000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.47|
|Biohydrogenation_with_weights|0.48|
|COVID_m1|3.73|
|COVID_m1_with_weights|1.57|
|ChemicalReactionNetwork|1.32|
|ChemicalReactionNetwork_with_weights|1.11|
|Cholera|203.83|
|Cholera_with_weights|219.35|
|DAISY_ex3|0.10|
|DAISY_ex3_with_weights|0.05|
|DAISY_mamil3|0.10|
|DAISY_mamil3_with_weights|0.09|
|DAISY_mamil4|3.34|
|DAISY_mamil4_with_weights|10.71|
|Goodwin| - |
|Goodwin_with_weights|711.06|
|HIV|0.33|
|HIV2|14.94|
|HIV2_with_weights|22.68|
|HIV_with_weights|0.21|
|LV|0.02|
|LV_with_weights|0.02|
|Lipolysis|0.12|
|Lipolysis_with_weights|0.10|
|NFkB|2049.79|
|NFkB_with_weights|971.21|
|OralGlucose|0.06|
|OralGlucose_with_weights|0.02|
|Pharm| - |
|Pharm_with_weights|2339.03|
|SEIR|0.18|
|SEIR2|0.11|
|SEIR2_with_weights|0.16|
|SEIRP|330.92|
|SEIRP_with_weights|84.65|
|SEIR_with_weights|0.08|
|SIRSForced|3.51|
|SIRSForced_with_weights|1.41|
|SIR_R0|0.01|
|SIR_R0_with_weights|0.01|
|SlowFast|0.04|
|SlowFast_with_weights|0.02|
|Treatment|0.12|
|Treatment_with_weights|0.13|
|Tumor|0.48|
|Tumor_with_weights|0.42|

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
