## Benchmark results

2024-01-12T08:15:25.284

Benchmarked backend: maple

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 3600 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.01|
|Biohydrogenation_with_weights|0.01|
|COVID_m1|0.72|
|COVID_m1_with_weights|0.14|
|ChemicalReactionNetwork|0.27|
|ChemicalReactionNetwork_with_weights|0.23|
|Cholera|42.16|
|Cholera_with_weights|46.12|
|DAISY_ex3|0.01|
|DAISY_ex3_with_weights|0.01|
|DAISY_mamil3|0.01|
|DAISY_mamil3_with_weights|0.01|
|DAISY_mamil4|0.42|
|DAISY_mamil4_with_weights|0.43|
|Goodwin|822.12|
|Goodwin_with_weights|146.60|
|HIV|0.02|
|HIV2|3.01|
|HIV2_with_weights|1.38|
|HIV_with_weights|0.01|
|LV|0.01|
|LV_with_weights|0.01|
|Lipolysis|0.01|
|Lipolysis_with_weights|0.01|
|NFkB|314.56|
|NFkB_with_weights|163.44|
|OralGlucose|0.01|
|OralGlucose_with_weights|0.01|
|Pharm| - |
|Pharm_with_weights|438.33|
|SEIR|0.01|
|SEIR2|0.03|
|SEIR2_with_weights|0.03|
|SEIRP|85.14|
|SEIRP_with_weights|19.86|
|SEIR_with_weights|0.01|
|SIRSForced|0.74|
|SIRSForced_with_weights|0.13|
|SIR_R0|0.01|
|SIR_R0_with_weights|0.01|
|SlowFast|0.01|
|SlowFast_with_weights|0.01|
|Treatment|0.01|
|Treatment_with_weights|0.02|
|Tumor|0.05|
|Tumor_with_weights|0.03|

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
