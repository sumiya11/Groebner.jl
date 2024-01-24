## Benchmark results

2024-01-23T15:40:45.178

Benchmarked backend: maple

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 300 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.01|
|Biohydrogenation_with_weights|0.01|
|COVID_m1|0.78|
|COVID_m1_with_weights|0.13|
|ChemicalReactionNetwork|0.27|
|ChemicalReactionNetwork_with_weights|0.22|
|Cholera|57.95|
|Cholera_with_weights|63.93|
|DAISY_ex3|0.01|
|DAISY_ex3_with_weights|0.02|
|DAISY_mamil3|0.01|
|DAISY_mamil3_with_weights|0.01|
|DAISY_mamil4|0.43|
|DAISY_mamil4_with_weights|0.45|
|Goodwin| - |
|Goodwin_with_weights| - |
|HIV|0.02|
|HIV2|3.54|
|HIV2_with_weights|1.62|
|HIV_with_weights|0.02|
|LV|0.01|
|LV_with_weights|0.01|
|Lipolysis|0.01|
|Lipolysis_with_weights|0.01|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose|0.01|
|OralGlucose_with_weights|0.01|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.01|
|SEIR2|0.06|
|SEIR2_with_weights|0.07|
|SEIRP|190.97|
|SEIRP_with_weights|28.67|
|SEIR_with_weights|0.01|
|SIRSForced|0.89|
|SIRSForced_with_weights|0.12|
|SIR_R0|0.01|
|SIR_R0_with_weights|0.01|
|SlowFast|0.01|
|SlowFast_with_weights|0.02|
|Treatment|0.02|
|Treatment_with_weights|0.02|
|Tumor|0.06|
|Tumor_with_weights|0.04|

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
