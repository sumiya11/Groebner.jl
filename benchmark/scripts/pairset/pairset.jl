using LinuxPerf, JLD2, BenchmarkTools

function clear_data()
    for (root, dirs, files) in walkdir((@__DIR__))
        for file in files
            if endswith((@__DIR__)*"/"*file, "jld2")
                rm((@__DIR__)*"/"*file)
            end
        end
        break
    end 
end

function load_data()
    data = []
    for (root, dirs, files) in walkdir((@__DIR__))
        @info "" files
        for file in files
            !endswith(file, "jld2") && continue
            a = load((@__DIR__)*"/"*file)
            push!(data, a)
        end
        break
    end 
    data
end

function describe_data(data)
    println("$(length(data)) entries")
    println("#\tb:filled\tps:load\t\tupdht:load\tupdht:sz\tpairs")
    println("-------------------------------------------------------------------------")
    res = []
    for (j, x) in enumerate(data)
        ht, update_ht, basis, pairset, i = x["ht"], x["update_ht"], x["basis"], x["pairset"], x["i"]
        # println("$j\t$(basis.nfilled)\t\t$(pairset.load)\t\t$(update_ht.load)/$(update_ht.size)")
        push!(res, [basis.nfilled, pairset.load, update_ht.load, update_ht.size, i])
    end
    println("avg\t", join(string.(round.(Int, mean(res, dims=1)[1])), "\t\t"))
    println("tot\t", join(string.(sum(res, dims=1)[1]), "\t\t"))
end

function multiple_data(data, n=1)
    newdata = deepcopy(data)
    for i in 1:n
        append!(newdata, deepcopy(data))
    end
    newdata
end

clear_data()

# data = load_data();
# describe_data(data)

_data = deepcopy(data);
@time begin
    for x in data
        ht, update_ht, basis, pairset, i = x["ht"], x["update_ht"], x["basis"], x["pairset"], x["i"]
        Groebner.pairset_update!(pairset, basis, ht, update_ht, i)
    end
end
#  0.005203 seconds (187 allocations: 5.844 KiB)
#  0.005697 seconds (187 allocations: 5.844 KiB)
#  0.005263 seconds (187 allocations: 5.844 KiB)

# _data = multiple_data(data);
_data = deepcopy(data);
@pstats "(cpu-cycles,task-clock),(instructions,branch-instructions,branch-misses), (L1-dcache-load-misses, L1-dcache-loads, cache-misses, cache-references)" begin
    for x in _data
        ht, update_ht, basis, pairset, i = x["ht"], x["update_ht"], x["basis"], x["pairset"], x["i"]
        Groebner.pairset_update!(pairset, basis, ht, update_ht, i)
    end
end

##############################

function gather_pairs(x::Vector{T}, pairs) where {T}
    s = T(0)
    @inbounds for p in pairs
        i, j = p
        s += x[i] + x[j]
    end
    s
end

n, m = 1000, 1000*1000
x = rand(UInt64, n)

pairs_i8 = [(rand(1:n), rand(1:n)) for i in 1:m]
pairs_i4 = [(Int32.(rand(1:n)), Int32.(rand(1:n))) for i in 1:m]
pairs_i2 = [(Int16.(rand(1:n)), Int16.(rand(1:n))) for i in 1:m]

@code_llvm debuginfo=:none gather_pairs(x, pairs_i8)
@code_native debuginfo=:none gather_pairs(x, pairs_i8)
@btime gather_pairs($x, pairs_i8) setup=begin
    pairs_i8 = [(Int64.(rand(1:n)), Int64.(rand(1:n))) for i in 1:m]
end

@code_llvm debuginfo=:none gather_pairs(x, pairs_i4)
@code_native debuginfo=:none gather_pairs(x, pairs_i4)
@btime gather_pairs($x, pairs_i4) setup=begin
    pairs_i4 = [(Int32.(rand(1:n)), Int32.(rand(1:n))) for i in 1:m]
end

@code_llvm debuginfo=:none gather_pairs(x, pairs_i2)
@code_native debuginfo=:none gather_pairs(x, pairs_i2)
@btime gather_pairs($x, pairs_i2) setup=begin
    pairs_i2 = [(Int16.(rand(1:n)), Int16.(rand(1:n))) for i in 1:m]
end

@pstats "(cpu-cycles,task-clock),(instructions,branch-instructions,branch-misses), (L1-dcache-load-misses, L1-dcache-loads, cache-misses, cache-references), (alignment-faults,page-faults,minor-faults)" begin
    for _ in 1:100
        gather_pairs(x, pairs_i8)
    end
end

@pstats "(cpu-cycles,task-clock),(instructions,branch-instructions,branch-misses), (L1-dcache-load-misses, L1-dcache-loads, cache-misses, cache-references), (alignment-faults,page-faults,minor-faults)" begin
    for _ in 1:100
        gather_pairs(x, pairs_i4)
    end
end

@pstats "(cpu-cycles,task-clock),(instructions,branch-instructions,branch-misses), (L1-dcache-load-misses, L1-dcache-loads, cache-misses, cache-references), (alignment-faults,page-faults,minor-faults)" begin
    for _ in 1:100
        gather_pairs(x, pairs_i2)
    end
end
