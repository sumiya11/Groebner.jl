using Base.Threads, BenchmarkTools

const COUNT = Ref{Int}(5 * 10^4)
const P = Ref{Float64}(0.05)

struct Monom
    value::Int
    visited::Bool
end

mutable struct Pool
    monoms::Vector{Monom}
    load::Int
end

function add_to_pool!(pool, value)
    @assert pool.load + 1 <= length(pool.monoms)
    pool.monoms[pool.load + 1] = Monom(value, false)
    pool.load += 1
    true
end

function pool_contains(pool, value)
    for i in 1:(pool.load)
        if pool.monoms[i].value == value
            return true
        end
    end
    false
end

function find_reducer!(basis, pool, i, max_elems)
    for j in 1:length(basis)
        if rand() < P[]
            salt = 0
            for i in 1:10
                salt = (salt + rand(1:(10^5))) % COUNT[]
                salt += floor(Int, sin(rand()) * cos(sin(rand())))
                salt += floor(Int, exp(rand()))
                salt += sum(i -> i^2, 1:rand(1:200))
            end
            sugar = pool.monoms[rand(1:(pool.load))].value
            value = powermod(salt, sugar, COUNT[])
            value = value % max_elems
            # @info "Generated $value"
            if !pool_contains(pool, value)
                add_to_pool!(pool, value)
            end
        end
    end
end

function symbolic_preprocessing!(basis, pool, max_elems=length(pool.monoms))
    i = 1
    calls = 0
    while i <= pool.load
        if pool.monoms[i].visited
            i += 1
            continue
        end

        pool.monoms[i] = Monom(pool.monoms[i].value, true)
        find_reducer!(basis, pool, i, max_elems)
        i += 1
        calls += 1
    end
    # @info "Final load = $(pool.load)"
    calls
end

##########################################################
##########################################################

function symbolic_preprocessing_threaded_1!(basis, pool)
    nchunks = nthreads()
    n_per_chunk = div(pool.load, nchunks)
    chunks = [
        pool.monoms[i:min(i + n_per_chunk - 1, pool.load)] for
        i in 1:n_per_chunk:(pool.load)
    ]
    el_per_chunk = map(length, chunks)
    # @info "Chunks: count = $nchunks, sizes = $(map(length, chunks))"
    pools = map(
        eldata -> Pool(resize!(eldata[2], div(COUNT[], 4)), eldata[1]),
        zip(el_per_chunk, chunks)
    )

    calls = Atomic{Int}()
    @threads for subpool in pools
        atomic_add!(calls, symbolic_preprocessing!(basis, subpool))
    end
    calls
end

##########################################################
##########################################################

@info "" nthreads()

n = max(div(COUNT[], 100), 1)
begin
    monoms = resize!(map(i -> Monom(i, 0), 1:n), COUNT[])
    pool1 = Pool(deepcopy(monoms), n)
    pool2 = Pool(deepcopy(monoms), n)
    basis1 = map(_ -> Monom(rand(1:COUNT[]), 0), 1:n)
    basis2 = map(_ -> Monom(rand(1:COUNT[]), 0), 1:n)

    symbolic_preprocessing!(basis1, pool1)
    symbolic_preprocessing_threaded_1!(basis2, pool2)
end

@btime symbolic_preprocessing!(basis, pool) setup = ((basis, pool) = begin
    monoms = resize!(map(i -> Monom(i, 0), 1:($n)), COUNT[])
    pool = Pool(deepcopy(monoms), $n)
    basis = map(_ -> Monom(rand(1:COUNT[]), 0), 1:($n))
    basis, pool
end)

@btime symbolic_preprocessing_threaded_1!(basis, pool) setup = ((basis, pool) = begin
    monoms = resize!(map(i -> Monom(i, 0), 1:($n)), COUNT[])
    pool = Pool(deepcopy(monoms), $n)
    basis = map(_ -> Monom(rand(1:COUNT[]), 0), 1:($n))
    basis, pool
end)
