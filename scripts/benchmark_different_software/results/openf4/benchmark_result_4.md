## Benchmark results

2024-01-12T00:10:51.376

Benchmarked backend: openf4

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.06|
|Biohydrogenation_with_weights| - |
|ChemicalReactionNetwork|2.30|
|ChemicalReactionNetwork_with_weights|14.89|
|Cholera| - |
|Cholera_with_weights| - |
|DAISY_ex3|0.03|
|DAISY_ex3_with_weights|33.04|
|DAISY_mamil3|0.03|
|DAISY_mamil3_with_weights|1.39|
|DAISY_mamil4|2.03|
|DAISY_mamil4_with_weights|43.27|
|HIV|1.13|
|HIV2|95.36|
|HIV2_with_weights| - |
|HIV_with_weights|945.47|
|LV|0.02|
|LV_with_weights|0.31|
|Lipolysis|0.02|
|Lipolysis_with_weights|0.07|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose|0.02|
|OralGlucose_with_weights|0.91|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.02|
|SEIR2|0.08|
|SEIR2_with_weights| - |
|SEIR_with_weights|2.95|
|SIRSForced|7.44|
|SIRSForced_with_weights| - |
|SIR_R0|0.02|
|SIR_R0_with_weights|0.02|
|SlowFast|0.02|
|SlowFast_with_weights|0.03|
|Treatment|0.07|
|Treatment_with_weights|12.35|
|Tumor|1.07|
|Tumor_with_weights| - |

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
