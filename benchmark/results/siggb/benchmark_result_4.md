## Benchmark results

2024-01-13T16:54:39.545

Benchmarked backend: siggb

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation| - |
|Biohydrogenation_with_weights| - |
|COVID_m1| - |
|COVID_m1_with_weights| - |
|ChemicalReactionNetwork| - |
|ChemicalReactionNetwork_with_weights| - |
|Cholera| - |
|Cholera_with_weights| - |
|DAISY_ex3| - |
|DAISY_ex3_with_weights| - |
|DAISY_mamil3|13.06|
|DAISY_mamil3_with_weights|11.90|
|DAISY_mamil4| - |
|DAISY_mamil4_with_weights| - |
|Goodwin| - |
|Goodwin_with_weights| - |
|HIV| - |
|HIV2| - |
|HIV2_with_weights| - |
|HIV_with_weights| - |
|LV|0.01|
|LV_with_weights|0.01|
|Lipolysis|0.23|
|Lipolysis_with_weights|0.14|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose| - |
|OralGlucose_with_weights| - |
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.44|
|SEIR2| - |
|SEIR2_with_weights| - |
|SEIRP| - |
|SEIRP_with_weights| - |
|SEIR_with_weights|0.19|
|SIRSForced| - |
|SIRSForced_with_weights| - |
|SIR_R0|0.06|
|SIR_R0_with_weights|0.07|
|SlowFast|0.01|
|SlowFast_with_weights|0.02|
|Treatment| - |
|Treatment_with_weights| - |
|Tumor| - |
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
