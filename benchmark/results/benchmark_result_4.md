## Benchmark results

2024-01-12T00:11:00.547

Benchmarked backends: Any["groebner", "maple", "msolve", "openf4", "singular"]

Benchmark suite: SIAN modulo 2^30 + 3

- Workers: 8
- Timeout: 1800 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|groebner|maple|msolve|openf4|singular|
|:----|---|---|---|---|---|
|Biohydrogenation|0.00|0.01|0.02|0.06|0.00|
|Biohydrogenation_with_weights|0.01|0.01|0.02| - |0.00|
|ChemicalReactionNetwork|0.34|0.29|0.61|2.30|51.54|
|ChemicalReactionNetwork_with_weights|0.61|0.20| - |14.89|42.96|
|Cholera|120.25|61.97|137.09| - | - |
|Cholera_with_weights|125.47|73.65|113.41| - | - |
|DAISY_ex3|0.00|0.02|0.02|0.03|0.00|
|DAISY_ex3_with_weights|0.00|0.01|15.56|33.04|0.00|
|DAISY_mamil3|0.00|0.02|0.02|0.03|0.00|
|DAISY_mamil3_with_weights|0.00|0.02|2.10|1.39|0.00|
|DAISY_mamil4|0.32|0.47|0.67|2.03|11.55|
|DAISY_mamil4_with_weights|0.32|0.47| - |43.27| - |
|HIV|0.01|0.02|0.03|1.13|0.01|
|HIV2|3.80|3.74|13.10|95.36|206.13|
|HIV2_with_weights|2.89|1.56|3.70| - |116.11|
|HIV_with_weights|0.01|0.01|0.02|945.47|0.01|
|LV|0.00|0.01|0.02|0.02|0.00|
|LV_with_weights|0.00|0.01|0.02|0.31|0.00|
|Lipolysis|0.00|0.01|0.01|0.02|0.00|
|Lipolysis_with_weights|0.00|0.01|0.02|0.07|0.00|
|NFkB| - |1137.15|1578.72| - | - |
|NFkB_with_weights|660.92|442.89| - | - | - |
|OralGlucose|0.00|0.02|0.02|0.02|0.00|
|OralGlucose_with_weights|0.00|0.01|0.02|0.91|0.00|
|Pharm| - | - | - | - | - |
|Pharm_with_weights| - |1612.93| - | - | - |
|SEIR|0.00|0.01|0.02|0.02|0.00|
|SEIR2|0.01|0.05|0.03|0.08|0.03|
|SEIR2_with_weights|0.08|0.06| - | - |0.03|
|SEIR_with_weights|0.01|0.02|0.08|2.95|0.00|
|SIRSForced|0.83|0.83|2.31|7.44|21.27|
|SIRSForced_with_weights|0.27|0.14|0.33| - |3.77|
|SIR_R0|0.00|0.01|0.02|0.02|0.00|
|SIR_R0_with_weights|0.00|0.01|0.02|0.02|0.00|
|SlowFast|0.00|0.01|0.02|0.02|0.00|
|SlowFast_with_weights|0.00|0.01|0.02|0.03|0.00|
|Treatment|0.01|0.02|0.03|0.07|0.01|
|Treatment_with_weights|0.07|0.02|0.03|12.35|0.01|
|Tumor|0.05|0.06|0.12|1.07|10.01|
|Tumor_with_weights|0.25|0.04|0.07| - | - |

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
