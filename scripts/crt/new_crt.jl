# Table of big integers must be initialized.
@timeit _TIMER function crt_vec_full!(
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    tables_ff::Vector{Vector{Vector{T}}},
    moduli::Vector{U},
    mask::Vector{BitVector}
) where {T <: Integer, U <: Integer}
    @invariant length(tables_ff) == length(moduli)
    @invariant length(table_zz) == length(mask)

    NN = Threads.threadpoolsize() + 1
    buf = [BigInt() for _ in 1:NN]
    n1 = [BigInt() for _ in 1:NN]
    n2 = [BigInt() for _ in 1:NN]

    mults = [Vector{BigInt}(undef, length(moduli)) for _ in 1:NN]
    for i in 1:length(mults[1])
        for t in 1:NN
            mults[t][i] = BigInt(0)
        end
    end
    # buf, n1, n2 = BigInt(), BigInt(), BigInt()
    # mults = Vector{BigInt}(undef, length(moduli))
    # for i in 1:length(mults)
    #     mults[i] = BigInt(0)
    # end

    crt_precompute!(modulo, n1[1], n2[1], mults[1], map(UInt64, moduli))
    for t in 2:NN
        Base.GMP.MPZ.set!(n1[t], n1[1])
        Base.GMP.MPZ.set!(n2[t], n2[1])
        for i in 1:length(mults[1])
            Base.GMP.MPZ.set!(mults[t][i], mults[1][i])
        end
    end

    # for t in 1:NN
    #     crt_precompute!(modulo, n1[t], n2[t], mults[t], map(UInt64, moduli))
    # end
    # crt_precompute!(modulo, n1, n2, mults, map(UInt64, moduli))

    @info "[crt] number of threads: $(NN)"

    # ind = Random.shuffle(1:length(table_zz)) # collect(1:length(table_zz)) 
    rems = [Vector{UInt64}(undef, length(moduli)) for _ in 1:NN]
    # rems = Vector{UInt64}(undef, length(moduli))
    @time for i in 1:length(table_zz)
        t = Threads.threadid()
        # Core.print("[thread $t] processing row $i / $(length(table_zz))\n")
        for j in 1:length(table_zz[i])
            mask[i][j] && continue

            for k in 1:length(moduli)
                @invariant 0 <= tables_ff[k][i][j] < moduli[k]
                rems[t][k] = UInt64(tables_ff[k][i][j])
            end

            crt!(modulo, buf[t], n1[t], n2[t], rems[t], mults[t])

            Base.GMP.MPZ.set!(table_zz[i][j], buf[t])
        end
    end
end
