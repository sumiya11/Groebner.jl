## Benchmark results

2024-02-02T17:56:45.998

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 10
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|Biohydrogenation|0.00|0.18|0.02|
|Biohydrogenation_with_weights|0.00|0.17|0.01|
|COVID_m1|0.50|0.78|1.72|
|COVID_m1_with_weights|0.12|0.31|0.30|
|ChemicalReactionNetwork|0.23|0.41|0.55|
|ChemicalReactionNetwork_with_weights|0.24|0.34| - |
|Cholera|43.00|46.54|111.14|
|Cholera_with_weights|45.60|50.90|102.71|
|DAISY_ex3|0.00|0.16|0.02|
|DAISY_ex3_with_weights|0.00|0.17|12.72|
|DAISY_mamil3|0.02|0.18|0.02|
|DAISY_mamil3_with_weights|0.00|0.16|1.83|
|DAISY_mamil4|0.22|0.55|0.61|
|DAISY_mamil4_with_weights|0.23|0.61| - |
|Goodwin| - | - | - |
|Goodwin_with_weights|152.54|255.78|458.67|
|HIV|0.01|0.17|0.02|
|HIV2|2.27|3.93|11.74|
|HIV2_with_weights|1.22|1.77|3.27|
|HIV_with_weights|0.00|0.20|0.02|
|LV|0.00|0.15|0.01|
|LV_with_weights|0.02|0.17|0.01|
|Lipolysis|0.00|0.20|0.02|
|Lipolysis_with_weights|0.00|0.21|0.02|
|NFkB|690.76|1057.66|1414.47|
|NFkB_with_weights|171.85|361.50| - |
|OralGlucose|0.00|0.21|0.02|
|OralGlucose_with_weights|0.00|0.20|0.02|
|Pharm| - | - | - |
|Pharm_with_weights|742.05|1495.19| - |
|SEIR|0.03|0.20|0.02|
|SEIR2|0.00|0.22|0.03|
|SEIR2_with_weights|0.02|0.27| - |
|SEIRP|76.15|130.81|160.54|
|SEIRP_with_weights|25.77|35.40|43.71|
|SEIR_with_weights|0.00|0.25|0.07|
|SIRSForced|0.63|1.29|2.36|
|SIRSForced_with_weights|0.13|0.33|0.35|
|SIR_R0|0.00|0.22|0.02|
|SIR_R0_with_weights|0.00|0.22|0.02|
|SlowFast|0.00|0.18|0.02|
|SlowFast_with_weights|0.00|0.21|0.02|
|Treatment|0.00|0.21|0.03|
|Treatment_with_weights|0.01|0.34|0.03|
|Tumor|0.04|0.26|0.11|
|Tumor_with_weights|0.06|0.22|0.08|

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
