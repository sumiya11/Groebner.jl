## Benchmark results

2024-01-24T05:02:24.970

Benchmarked backends: Any["groebner", "maple", "msolve"]

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 16
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|
|:----|---|---|---|
|Biohydrogenation|0.00|0.01|0.02|
|Biohydrogenation_with_weights|0.00|0.01|0.02|
|COVID_m1|0.52|0.73|1.72|
|COVID_m1_with_weights|0.13|0.13|0.31|
|ChemicalReactionNetwork|0.22|0.27|0.54|
|ChemicalReactionNetwork_with_weights|0.25|0.21| - |
|Cholera|54.74|57.19|132.39|
|Cholera_with_weights|60.09|61.35|114.62|
|DAISY_ex3|0.00|0.01|0.01|
|DAISY_ex3_with_weights|0.00|0.01|14.29|
|DAISY_mamil3|0.02|0.01|0.02|
|DAISY_mamil3_with_weights|0.00|0.01|1.83|
|DAISY_mamil4|0.25|0.44|0.68|
|DAISY_mamil4_with_weights|0.31|0.45| - |
|Goodwin| - | - | - |
|Goodwin_with_weights|428.86|354.33|599.37|
|HIV|0.01|0.02|0.03|
|HIV2|2.31|3.56|12.02|
|HIV2_with_weights|1.32|1.62|4.18|
|HIV_with_weights|0.00|0.01|0.02|
|LV|0.00|0.01|0.01|
|LV_with_weights|0.00|0.01|0.01|
|Lipolysis|0.00|0.01|0.02|
|Lipolysis_with_weights|0.00|0.02|0.02|
|NFkB| - | - |1763.82|
|NFkB_with_weights|505.92|453.74| - |
|OralGlucose|0.02|0.02|0.02|
|OralGlucose_with_weights|0.00|0.02|0.02|
|Pharm| - | - | - |
|Pharm_with_weights| - | - | - |
|SEIR|0.05|0.01|0.02|
|SEIR2|0.00|0.07|0.03|
|SEIR2_with_weights|0.02|0.07| - |
|SEIRP|153.94|190.14|212.19|
|SEIRP_with_weights|28.59|28.39|67.67|
|SEIR_with_weights|0.00|0.02|0.07|
|SIRSForced|0.72|0.90|1.95|
|SIRSForced_with_weights|0.17|0.12|0.32|
|SIR_R0|0.00|0.01|0.02|
|SIR_R0_with_weights|0.00|0.01|0.02|
|SlowFast|0.00|0.01|0.02|
|SlowFast_with_weights|0.00|0.02|0.02|
|Treatment|0.01|0.02|0.03|
|Treatment_with_weights|0.01|0.02|0.03|
|Tumor|0.08|0.06|0.13|
|Tumor_with_weights|0.09|0.04|0.07|

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
