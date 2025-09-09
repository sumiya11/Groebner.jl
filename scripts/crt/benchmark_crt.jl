using Base.Threads

function crt!(M, buf, n1, n2, ai, ci)
    Base.GMP.MPZ.set_ui!(n1, UInt(0))
    for i in 1:length(ai)
        Base.GMP.MPZ.mul_ui!(n2, ci[i], ai[i])
        Base.GMP.MPZ.add!(n1, n2)
    end
    Base.GMP.MPZ.set!(buf, n1)
    Base.GMP.MPZ.fdiv_r!(buf, M)
end

function crt_vec_full!(table_zz, modulo, tables_ff, mults, moduli)
    NN = Threads.maxthreadid()
    buf = [BigInt() for _ in 1:NN]
    n1 = [BigInt() for _ in 1:NN]
    n2 = [BigInt() for _ in 1:NN]
    rems = [Vector{UInt64}(undef, length(moduli)) for _ in 1:NN]
    Threads.@threads :static for i in 1:length(table_zz)
        t = Threads.threadid()
        for j in 1:length(table_zz[i])
            for k in 1:length(moduli)
                rems[t][k] = UInt64(tables_ff[k][i][j])
            end
            crt!(modulo, buf[t], n1[t], n2[t], rems[t], mults)
            Base.GMP.MPZ.set!(table_zz[i][j], buf[t])
        end
    end
    @assert NN == Threads.maxthreadid()
end

# Preallocate data
M, P = 512, 101
table_zz = [[BigInt() for _ in 1:M] for i in 1:M]
moduli = UInt64.(prevprimes(2^31-1, P))
modulo = prod(BigInt, moduli)
mults = [rand(1:modulo-1) for _ in 1:length(moduli)]
tables_ff = [ [[mod(rand(UInt64), moduli[j]) for _ in 1:M] for i in 1:M] for j in 1:P]

@time crt_vec_full!(table_zz,modulo,tables_ff,mults,moduli)
@time crt_vec_full!(table_zz,modulo,tables_ff,mults,moduli)
