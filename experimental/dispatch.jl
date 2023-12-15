using AbstractAlgebra, AllocCheck, BenchmarkTools

R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)

s = [x * y^2 + z, x^2 * z + y]
S = Groebner.katsuran(9, ground=GF(2^31 - 1), ordering=:degrevlex)

####################################
############  DEFAULT  #############
####################################

Groebner.logging_enabled() = true
r = check_allocs(Groebner.groebner, (typeof(s),))
@info length(r)  # 132
@benchmark Groebner.groebner($s)
#=
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  134.600 μs …   8.695 ms  ┊ GC (min … max): 0.00% … 94.07%     
 Time  (median):     148.800 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   162.509 μs ± 269.036 μs  ┊ GC (mean ± σ):  5.49% ±  3.27%

  ▂▇█▇▆▅▇▇▆▅▄▃▃▂▁ ▁  ▁                                          ▂
  █████████████████████▇███▇▆▇▇▇▇▆▇▇▆▇▇▅▆▆▆▆▅▆▅▆▅▆▅▅▅▆▅▅▅▅▆▄▅▅▄ █
  135 μs        Histogram: log(frequency) by time        276 μs <

 Memory estimate: 55.59 KiB, allocs estimate: 581.
=#

@benchmark Groebner.groebner($S)
#=
BenchmarkTools.Trial: 23 samples with 1 evaluation.
 Range (min … max):  197.985 ms … 304.716 ms  ┊ GC (min … max): 2.03% … 3.13%
 Time  (median):     213.912 ms               ┊ GC (median):    1.88%
 Time  (mean ± σ):   221.748 ms ±  24.122 ms  ┊ GC (mean ± σ):  2.16% ± 0.66%      

     ▃▃█   ▃ █   
  ▇▁▁███▇▇▇█▁█▁▁▇▁▁▁▇▁▇▇▁▁▁▁▁▁▁▁▁▁▇▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  198 ms           Histogram: frequency by time          305 ms <

 Memory estimate: 54.83 MiB, allocs estimate: 35010.
=#

####################################
############  NO LOGS  #############
####################################

Groebner.logging_enabled() = false
r = check_allocs(Groebner.groebner, (typeof(s),))
@info length(r)  # 66
@benchmark Groebner.groebner($s)
#=
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  35.500 μs …   7.040 ms  ┊ GC (min … max):  0.00% … 95.62%
 Time  (median):     39.500 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   49.397 μs ± 201.520 μs  ┊ GC (mean ± σ):  11.76% ±  2.87%     

   ▃█▃           
  ▂████▆▄▃▂▂▂▃▃▃▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  35.5 μs         Histogram: frequency by time         96.8 μs <

 Memory estimate: 42.92 KiB, allocs estimate: 308.
=#

@benchmark Groebner.groebner($S)
#=
BenchmarkTools.Trial: 22 samples with 1 evaluation.
 Range (min … max):  194.198 ms … 391.634 ms  ┊ GC (min … max): 3.33% … 2.89%
 Time  (median):     226.074 ms               ┊ GC (median):    1.91%
 Time  (mean ± σ):   234.905 ms ±  44.547 ms  ┊ GC (mean ± σ):  2.12% ± 0.67%      

  ▃  ▃ ▃     ▃ █       ▃
  █▇▇█▇█▁▇▇▁▁█▇█▁▁▁▇▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  194 ms           Histogram: frequency by time          392 ms <

 Memory estimate: 54.80 MiB, allocs estimate: 34297.
=#

####################################
#######  FIXING DISPATCH  ##########
####################################

Groebner.logging_enabled() = false
r = check_allocs(Groebner.groebner, (typeof(s),))
@info length(r)  # 37
@benchmark Groebner.groebner($s)
#=
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  30.400 μs …  10.521 ms  ┊ GC (min … max):  0.00% … 98.74%
 Time  (median):     36.400 μs               ┊ GC (median):     0.00%    
 Time  (mean ± σ):   50.917 μs ± 272.105 μs  ┊ GC (mean ± σ):  15.76% ±  
2.95%

  ▆██▇▆▆▅▅▅▅▄▄▃▃▃▃▃▃▂▂▁▁▁▁  ▁▁▁▁▁▁                             ▂
  █████████████████████████████████▇▇▇▇▆█▇▇▇▆▇▆▅▅▆▆▆▆▄▅▅▆▅▅▅▄▅ █
  30.4 μs       Histogram: log(frequency) by time       119 μs <

 Memory estimate: 41.78 KiB, allocs estimate: 290.
=#

@benchmark Groebner.groebner($S)
#=
BenchmarkTools.Trial: 27 samples with 1 evaluation.
 Range (min … max):  178.942 ms … 231.366 ms  ┊ GC (min … max): 2.46% … 1.76%
 Time  (median):     186.839 ms               ┊ GC (median):    2.31%
 Time  (mean ± σ):   189.678 ms ±  11.078 ms  ┊ GC (mean ± σ):  2.58% ± 0.78%        

   ▃  ██   █  ▃
  ▇█▇▇██▁▇▇█▁▇█▇▇▁▁▁▁▇▁▁▇▁▇▁▁▁▇▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  179 ms           Histogram: frequency by time          231 ms <

 Memory estimate: 54.80 MiB, allocs estimate: 34279.
=#

@profview_allocs Groebner.groebner(s) sample_rate = 1.0
