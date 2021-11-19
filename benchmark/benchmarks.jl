
using BenchmarkTools

fs = GroebnerBases.katsura6(ground=GF(2^31-1))
@benchmark f4(fs)
