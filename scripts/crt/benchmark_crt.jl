using Revise
using Primes, Groebner

if false
M = 1024
P = 101
modulo = BigInt()
table_zz = [[BigInt() for _ in 1:i] for i in 1:M]
moduli = UInt64.(prevprimes(2^31-1, P))
tables_ff = [ [[mod(rand(UInt64), moduli[j]) for _ in 1:i] for i in 1:M] for j in 1:P]
mask = [falses(i) for i in 1:M]
end

@time Groebner.crt_vec_full!(
    table_zz,
    modulo,
    tables_ff,
    moduli,
    mask
)
