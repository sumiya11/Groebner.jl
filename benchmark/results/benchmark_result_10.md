## Benchmark results

2024-02-03T04:20:38.714

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: SIAN, the rationals

- Workers: 10
- Timeout: 3000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|Biohydrogenation|0.10|0.47|0.42|
|Biohydrogenation_with_weights|0.08|0.48|0.36|
|COVID_m1|2.87|3.73|10.85|
|COVID_m1_with_weights|0.87|1.57|8.36|
|ChemicalReactionNetwork|2.11|1.32|5.17|
|ChemicalReactionNetwork_with_weights|1.29|1.11|3.84|
|Cholera|1100.59|203.83|1841.30|
|Cholera_with_weights|1375.00|219.35|1639.68|
|DAISY_ex3|0.00|0.10|0.06|
|DAISY_ex3_with_weights|0.03|0.05|0.10|
|DAISY_mamil3|0.01|0.10|0.08|
|DAISY_mamil3_with_weights|0.03|0.09|0.11|
|DAISY_mamil4|1.20|3.34|2.72|
|DAISY_mamil4_with_weights|1.21|10.71|7.22|
|Goodwin| - | - | - |
|Goodwin_with_weights| - |711.06|2240.06|
|HIV|0.04|0.33|0.30|
|HIV2|16.27|14.94|67.04|
|HIV2_with_weights|9.14|22.68|41.33|
|HIV_with_weights|0.02|0.21|0.30|
|LV|0.00|0.02|0.04|
|LV_with_weights|0.00|0.02|0.04|
|Lipolysis|0.01|0.12|0.05|
|Lipolysis_with_weights|0.01|0.10|0.04|
|NFkB| - |2049.79| - |
|NFkB_with_weights| - |971.21| - |
|OralGlucose|0.00|0.06|0.04|
|OralGlucose_with_weights|0.00|0.02|0.07|
|Pharm| - | - | - |
|Pharm_with_weights| - |2339.03| - |
|SEIR|0.01|0.18|0.06|
|SEIR2|0.04|0.11|0.16|
|SEIR2_with_weights|0.10|0.16|0.44|
|SEIRP| - |330.92|2450.64|
|SEIRP_with_weights|216.98|84.65|313.82|
|SEIR_with_weights|0.01|0.08|0.07|
|SIRSForced|3.00|3.51|11.94|
|SIRSForced_with_weights|1.03|1.41|8.93|
|SIR_R0|0.00|0.01|0.04|
|SIR_R0_with_weights|0.00|0.01|0.05|
|SlowFast|0.00|0.04|0.04|
|SlowFast_with_weights|0.04|0.02|0.04|
|Treatment|0.02|0.12|0.14|
|Treatment_with_weights|0.03|0.13|0.18|
|Tumor|0.23|0.48|0.82|
|Tumor_with_weights|0.25|0.42|10.30|

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
