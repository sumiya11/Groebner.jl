## Benchmark results

2024-01-24T04:21:16.720

Benchmarked backend: groebner

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.00|
|Biohydrogenation_with_weights|0.00|
|COVID_m1|0.52|
|COVID_m1_with_weights|0.13|
|ChemicalReactionNetwork|0.22|
|ChemicalReactionNetwork_with_weights|0.25|
|Cholera|54.74|
|Cholera_with_weights|60.09|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.00|
|DAISY_mamil3|0.02|
|DAISY_mamil3_with_weights|0.00|
|DAISY_mamil4|0.25|
|DAISY_mamil4_with_weights|0.31|
|Goodwin| - |
|Goodwin_with_weights|428.86|
|HIV|0.01|
|HIV2|2.31|
|HIV2_with_weights|1.32|
|HIV_with_weights|0.00|
|LV|0.00|
|LV_with_weights|0.00|
|Lipolysis|0.00|
|Lipolysis_with_weights|0.00|
|NFkB| - |
|NFkB_with_weights|505.92|
|OralGlucose|0.02|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.05|
|SEIR2|0.00|
|SEIR2_with_weights|0.02|
|SEIRP|153.94|
|SEIRP_with_weights|28.59|
|SEIR_with_weights|0.00|
|SIRSForced|0.72|
|SIRSForced_with_weights|0.17|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.00|
|SlowFast|0.00|
|SlowFast_with_weights|0.00|
|Treatment|0.01|
|Treatment_with_weights|0.01|
|Tumor|0.08|
|Tumor_with_weights|0.09|

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
