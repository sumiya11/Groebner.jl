## Benchmark results

2024-02-03T06:59:32.798

Benchmarked backend: groebner

Benchmark suite: SIAN, the rationals

- Workers: 10
- Timeout: 7200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|Biohydrogenation|0.09|
|Biohydrogenation_with_weights|0.07|
|COVID_m1|2.62|
|COVID_m1_with_weights|0.78|
|ChemicalReactionNetwork|1.89|
|ChemicalReactionNetwork_with_weights|1.14|
|Cholera|1035.43|
|Cholera_with_weights|1353.44|
|DAISY_ex3|0.00|
|DAISY_ex3_with_weights|0.03|
|DAISY_mamil3|0.01|
|DAISY_mamil3_with_weights|0.01|
|DAISY_mamil4|1.10|
|DAISY_mamil4_with_weights|1.14|
|Goodwin| - |
|Goodwin_with_weights|1442.35|
|HIV|0.03|
|HIV2|14.69|
|HIV2_with_weights|8.43|
|HIV_with_weights|0.02|
|LV|0.00|
|LV_with_weights|0.00|
|Lipolysis|0.01|
|Lipolysis_with_weights|0.01|
|NFkB| - |
|NFkB_with_weights| - |
|OralGlucose|0.00|
|OralGlucose_with_weights|0.00|
|Pharm| - |
|Pharm_with_weights| - |
|SEIR|0.01|
|SEIR2|0.02|
|SEIR2_with_weights|0.09|
|SEIRP|1682.08|
|SEIRP_with_weights|215.62|
|SEIR_with_weights|0.01|
|SIRSForced|2.83|
|SIRSForced_with_weights|0.82|
|SIR_R0|0.00|
|SIR_R0_with_weights|0.00|
|SlowFast|0.00|
|SlowFast_with_weights|0.00|
|Treatment|0.02|
|Treatment_with_weights|0.03|
|Tumor|0.24|
|Tumor_with_weights|0.24|

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
