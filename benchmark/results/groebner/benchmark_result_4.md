## Benchmark results

2024-01-12T08:37:44.375

Benchmarked backend: groebner

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 1000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.00|
|Biohydrogenation_with_weights|0.01|
|COVID_m1|0.65|
|COVID_m1_with_weights|0.19|
|ChemicalReactionNetwork|0.38|
|ChemicalReactionNetwork_with_weights|0.55|
|Cholera|52.40|
|Cholera_with_weights|57.69|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.00|
|DAISY_mamil3|0.00|
|DAISY_mamil3_with_weights|0.00|
|DAISY_mamil4|0.32|
|DAISY_mamil4_with_weights|0.31|
|Goodwin| - |
|Goodwin_with_weights|223.73|
|HIV|0.01|
|HIV2|3.14|
|HIV2_with_weights|2.16|
|HIV_with_weights|0.03|
|LV|0.00|
|LV_with_weights|0.00|
|Lipolysis|0.00|
|Lipolysis_with_weights|0.00|
|NFkB| - |
|NFkB_with_weights|285.44|
|OralGlucose|0.00|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.00|
|SEIR2| - |
|SEIR2_with_weights|0.03|
|SEIRP| - |
|SEIRP_with_weights|34.90|
|SEIR_with_weights| - |
|SIRSForced|0.72|
|SIRSForced_with_weights|0.27|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.02|
|SlowFast| - |
|SlowFast_with_weights|0.00|
|Treatment| - |
|Treatment_with_weights|0.03|
|Tumor|0.05|
|Tumor_with_weights|0.18|

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
