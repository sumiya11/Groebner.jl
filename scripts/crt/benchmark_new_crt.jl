using Base.Threads: nthreads, @threads, @spawn
using Base.Iterators: partition
using Primes, BenchmarkTools

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
    # NN = Threads.maxthreadid()
    # buf = [BigInt() for _ in 1:NN]
    # n1 = [BigInt() for _ in 1:NN]
    # n2 = [BigInt() for _ in 1:NN]
    # rems = [Vector{UInt64}(undef, length(moduli)) for _ in 1:NN]
    tasks_per_thread = 1
    n_threads = nthreads()
    chunk_size = max(1, div(length(table_zz), (tasks_per_thread * n_threads)))
    data_chunks = partition(Random.shuffle(1:length(table_zz)), chunk_size)
    tasks = Task[]
    for chunk in data_chunks
        # @info "chunk" chunk
        task = @spawn begin
            # state = some_initial_value
            buf = BigInt()
            n1 = BigInt()
            n2 = BigInt()
            rems = Vector{UInt64}(undef, length(moduli))
            for i in chunk
                for j in 1:length(table_zz[i])
                    for k in 1:length(moduli)
                        rems[k] = UInt64(tables_ff[k][i][j])
                    end
                    crt!(modulo, buf, n1, n2, rems, mults)
                    Base.GMP.MPZ.set!(table_zz[i][j], buf)
                end
            end
            # return state
        end
        push!(tasks, task)
    end
    for task in tasks
        wait(task)
    end
    #=
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
    =#
end

# Preallocate data
M, P = 512, 101
table_zz = [[BigInt() for _ in 1:i] for i in 1:M]
moduli = UInt64.(prevprimes(2^31-1, P))
modulo = prod(BigInt, moduli)
mults = [rand(1:modulo-1) for _ in 1:length(moduli)]
tables_ff = [ [[mod(rand(UInt64), moduli[j]) for _ in 1:i] for i in 1:M] for j in 1:P]

@time crt_vec_full!(table_zz,modulo,tables_ff,mults,moduli)
@btime crt_vec_full!($table_zz,$modulo,$tables_ff,$mults,$moduli)
