## Benchmark results

2024-01-12T06:08:17.235

Benchmarked backend: groebner

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 5000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.00|
|Biohydrogenation_with_weights|0.01|
|ChemicalReactionNetwork|0.33|
|ChemicalReactionNetwork_with_weights|0.51|
|Cholera|65.40|
|Cholera_with_weights|69.96|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.00|
|DAISY_mamil3|0.00|
|DAISY_mamil3_with_weights|0.00|
|DAISY_mamil4|0.32|
|DAISY_mamil4_with_weights|0.34|
|HIV|0.01|
|HIV2|3.35|
|HIV2_with_weights|2.58|
|HIV_with_weights|0.01|
|LV|0.00|
|LV_with_weights|0.00|
|Lipolysis|0.00|
|Lipolysis_with_weights|0.00|
|NFkB|1629.72|
|NFkB_with_weights|559.35|
|OralGlucose|0.00|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights|1922.79|
|SEIR|0.00|
|SEIR2|0.01|
|SEIR2_with_weights|0.05|
|SEIR_with_weights|0.00|
|SIRSForced|0.76|
|SIRSForced_with_weights|0.20|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.00|
|SlowFast|0.00|
|SlowFast_with_weights|0.00|
|Treatment|0.01|
|Treatment_with_weights|0.03|
|Tumor|0.05|
|Tumor_with_weights|0.21|

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
