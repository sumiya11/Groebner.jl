## Benchmark results

<<<<<<< HEAD
2023-09-25T16:17:24.805

- Benchmarked backend: `groebner`
- Workers: 16
- Timeout: 300 s
=======
2023-09-25T18:07:16.254

- Benchmarked backend: `groebner`
- Workers: 4
- Timeout: 30 s
>>>>>>> 90f63b68f55062bea6bf46570f976f4d0d338637

**All timings in seconds.**

|Model|Total|
|-----|---|
<<<<<<< HEAD
|cyclic 7|0.10|
|cyclic 8|1.18|
|cyclic 9|96.56|
|dummy|0.00|
|eco 11|0.32|
|eco 12|1.84|
|eco 13|8.25|
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|1.75|
|katsura 10|0.67|
|katsura 11|4.50|
|katsura 12|31.70|
|noon 7|0.16|
|noon 8|1.52|
|noon 9|15.80|
|reimer 6|0.05|
|reimer 7|0.71|
|reimer 8|16.84|

*Benchmarking environment:*

* Total RAM (GiB): 2003
* Processor: AMD EPYC 7702 64-Core Processor                
=======
|cyclic 7|0.14|
|cyclic 8|1.80|
|cyclic 9| - |
|dummy|0.00|
|eco 11|0.69|
|eco 12| - |
|eco 13| - |
|henrion 5|0.00|
|henrion 6|0.05|
|henrion 7| - |
|katsura 10|1.52|
|katsura 11| - |
|katsura 12| - |
|noon 7|0.35|
|noon 8| - |
|noon 9| - |
|reimer 6|0.11|
|reimer 7|1.43|
|reimer 8| - |

*Benchmarking environment:*

* Total RAM (GiB): 7
* Processor: Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz
>>>>>>> 90f63b68f55062bea6bf46570f976f4d0d338637
* Julia version: 1.9.1

Versions of the dependencies:

* Combinatorics : 1.0.2
* MultivariatePolynomials : 0.5.2
* Primes : 0.5.4
* ExprTools : 0.1.10
* SIMD : 3.4.5
* AbstractAlgebra : 0.32.1
* SnoopPrecompile : 1.0.3
