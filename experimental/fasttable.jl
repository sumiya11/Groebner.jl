
using Plots
using Statistics
using Random

function degtype()
    UInt32
end

function monoms_of_degree_not_greater_than(N, degree)
    if degree == 0
        [zeros(degtype(), N)]
    else
        prev = monoms_of_degree_not_greater_than(N, degree - 1)
        one_at_i = (N, i) -> (x = zeros(degtype(), N); x[i] += 1; x)
        filter(
            x -> count(!iszero, x) < 4,
            union(prev, [[pk + one_at_i(N, i) for pk in prev] for i in 1:N]...)
        )
    end
end

function linear_hash_all_buckets(N, H, deg)
    h = rand(degtype(), N)
    h2 = rand(degtype())
    monoms = shuffle(monoms_of_degree_not_greater_than(N, deg))
    println("Hash vector:\n\t$h")
    println("Monoms: $(length(monoms))")
    println("Hashtable: $(length(H))")
    for m in monoms
        hm = sum(h .* m) + count(!iszero, m) * h2
        hm = hm % length(H) + 1
        if length(H[hm]) > 0
            # println("collision: $(H[hm]) and $m")
        end
        push!(H[hm], m)
    end
    println("Overhead >1: $(sum(length, filter(x -> length(x) > 1, H)))")
    H
end

function linear_hash_all_openaddr(N, H, deg)
    collided = 0
    colarray = []
    longestchain = 0
    h = rand(degtype(), N)
    h2 = rand(degtype())
    monoms = shuffle(monoms_of_degree_not_greater_than(N, deg))
    println("Hash vector:\n\t$h")
    println("Monoms: $(length(monoms))")
    println("Hashtable: $(length(H))")
    for m in monoms
        col = 0
        hm = sum(h .* m)
        hm = hm % length(H) + 1
        while length(H[hm]) > 0
            # @info "collided"
            collided += 1
            col += 1
            hm = hm % length(H) + 1
        end
        longestchain = max(longestchain, col)
        push!(colarray, col)
        push!(H[hm], m)
    end
    println("\nCollided: $collided")
    println("Longest chain: $longestchain")
    println("Chain average: $(maximum(colarray))")
    H, colarray
end

M = 2^15
H = [[] for i in 1:M]

N = 9
deg = 10

H = linear_hash_all_buckets(N, H, deg)

lengthes = Dict(i => 0 for i in 0:30)

for j in 1:length(H)
    lengthes[length(H[j])] += 1
end

# histogram(xxx)

plot(collect(1:length(H)), map(length, H[1:end]))
# plot(lengthes)
