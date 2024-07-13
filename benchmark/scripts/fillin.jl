using LinearAlgebra, SparseArrays
using Random

mutable struct MyMatrix
    supports::Vector{Vector{Any}}
    coeffs::Vector{Vector{Any}}
    m::Int
    n::Int
end

Base.size(mat::MyMatrix) = (mat.m, mat.n)

function matrix_from_dense(arr)
    inds = findall(!iszero, arr)
    supports = [findall(!iszero, arr[i, :]) for i in 1:size(arr, 1)]
    coeffs = [[arr[i, j] for j in supports[i]] for i in 1:size(arr, 1)]
    MyMatrix(supports, coeffs, size(arr)...)
end

function matrix_to_sparse(mat::MyMatrix)
    m, n = size(mat)
    arr = spzeros(Int, m, n)
    for i in 1:m
        for j in 1:length(mat.supports[i])
            row, col = i, mat.supports[i][j]
            arr[row, col] = 1
        end
    end
    sparse(arr)
end

function build_graph(mat)
    m, n = size(mat)
    graph = [Int[] for _ in 1:(m + n)]
    for i in 1:m
        for j in 1:length(mat.supports[i])
            row, col = i, mat.supports[i][j]
            push!(graph[row], m + col)
            push!(graph[m + col], row)
        end
    end
    graph
end

function dfs!(graph, src, components, color)
    if components[src] != 0
        return 0
    end
    queue = Int[src]
    while !isempty(queue)
        src = pop!(queue)
        components[src] = color
        for dst in graph[src]
            if components[dst] == 0
                push!(queue, dst)
            end
        end
    end
    return 1
end

function connected_components(graph)
    components = [0 for _ in 1:length(graph)]
    color = 1
    for i in 1:length(graph)
        color += dfs!(graph, i, components, color)
    end
    components
end

function decompose(mat::MyMatrix)
    graph = build_graph(mat)
    components = connected_components(graph)
    components
end

A = [
    1 0 1;
    0 1 0;
    0 0 1;
]
B = matrix_from_dense(A)
A = matrix_to_sparse(B)
cc = decompose(B)

begin
    n = 100
    A = sprand(Int, n, n, 0.2) + I
    for i in 1:n
        for j in (i+1):n
            A[j,i] = 0
        end
    end
    dropzeros!(A)
    @info "" A
    B = matrix_from_dense(A)
    cc = decompose(B)
    length(unique(cc)), map(i -> sum(cc .== i), 1:length(unique(cc)))
end

#=
Int32[2, 3, 4, 5]
UInt32[0x00000001, 0x40000002, 0x20000002, 0x20000001]

Int32[2, 3, 4, 5]
UInt32[0x00000001, 0x40000002, 0x20000002, 0x20000001]
=#

include("../../../Groebner.jl/src/Groebner.jl")
using AbstractAlgebra, JLD2, Random
c = Groebner.Examples.noonn(8, k=GF(2^30+3), internal_ordering=:degrevlex)

Groebner.__SAVE[] = true
gb1 = Groebner.groebner(c, loglevel=-0, linalg=:deterministic);
Groebner.__SAVE[] = false
gb2 = Groebner.groebner(c, loglevel=-0, linalg=:experimental_2);
gb1 == gb2

for fn in readdir((@__DIR__), join=true)
    if !startswith(last(split(fn, "/")), "matrix")
        continue
    end
    begin
        mat = load(fn)["matrix"]
        (n, _) = Groebner.matrix_block_sizes(mat)
        rows = map(row -> filter(col -> col <= n, row), mat.upper_rows)
        # B = MyMatrix(rows, [[]], n, n)
        # A = matrix_to_sparse(B)
        # @time cc = decompose(B)
        # nc = length(unique(cc))
        # cnt = map(i -> floor(Int, sum(cc .== unique(cc)[i]) / 2), 1:length(unique(cc)))
        # display(A)
        @info "" n
    end
end

for fn in readdir((@__DIR__), join=true)
    if !startswith(last(split(fn, "/")), "matrix")
        continue
    end
    begin
        matrix = load(fn)["matrix"]
        basis = load(fn)["basis"]
        println(Groebner.matrix_block_sizes(matrix))
        arithmetic = Groebner.SpecializedArithmeticZp(UInt64, UInt32, 2^30 + 3)
        @time Groebner.linalg_deterministic_sparse!(
            matrix, 
            basis, 
            Groebner.LinearAlgebra(:a, :b), 
            arithmetic,
        )

        matrix = load(fn)["matrix"]
        basis = load(fn)["basis"]
        @time Groebner.linalg_randomized_sparse!(
            matrix, 
            basis, 
            Groebner.LinearAlgebra(:a, :b),
            arithmetic,
            Random.MersenneTwister(42)
        )

        matrix = load(fn)["matrix"]
        basis = load(fn)["basis"]
        @time Groebner.linalg_randomized_hashcolumns_sparse!(
            matrix, 
            basis, 
            Groebner.LinearAlgebra(:a, :b),
            arithmetic,
            Random.MersenneTwister(42)
        )
    end
end

if true
    for fn in readdir((@__DIR__), join=true)
        if !startswith(last(split(fn, "/")), "matrix")
            continue
        end
        rm(fn)
    end
end
