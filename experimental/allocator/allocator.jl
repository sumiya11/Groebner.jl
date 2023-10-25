using BenchmarkTools, Random
using StrideArrays

mutable struct BumpAllocator
    buf::Vector{UInt8}
    offset::UInt
end

function BumpAllocator(n::Integer)
    @debug "Allocating buffer of $n bytes"
    BumpAllocator(Vector{UInt8}(undef, n), 0)
end

# function alloc(b::BumpAllocator, ::Type{Array{T, N}}, d::NTuple{N, Int}) where {T, N}
#     # @debug "Allocating $(prod(d)) instances of $(T) ($(b.offset) / $(length(b.buf)) bytes filled)"
#     # b.offset + sizeof(T) * prod(d) > length(b.buf) && error("Alloc: Out of memory")
#     ptr = Base.unsafe_convert(Ptr{UInt8}, b.buf) + b.offset
#     b.offset += sizeof(T) * prod(d)
#     unsafe_wrap(Array, convert(Ptr{T}, ptr), d)
# end

function alloc(b::BumpAllocator, ::Type{Array{T, N}}, d::NTuple{N, Int}) where {T, N}
    ptr = pointer(b.buf) + b.offset
    sz = prod(d)
    b.offset += sz
    # b.offset > sizeof(b.buf) && oom_error(b)
    ptr_ = convert(Ptr{T}, ptr)
    PtrArray(ptr_, d)
end

function clear!(b::BumpAllocator)
    b.offset = 0
    return b
end

#########################
#########################

compute_1(a, b, c) = a^42 % (b * c + 1)
compute_2(a) = reduce(+, map(sum, a))

function do_stuff_1(m, n, len=1:100)
    s = 0
    for i in 1:m
        vecs = Vector{Vector{Int}}(undef, n)
        for j in 1:n
            vecs[j] = Vector{Int}(undef, rand(len))
            for k in 1:length(vecs[j])
                vecs[j][k] = compute_1(i, j, k)
            end
        end
        s = s + compute_2(vecs)
    end
    s
end

function do_stuff_2(m, n, len=1:100)
    s = 0
    allocator = BumpAllocator(n * maximum(len) * sizeof(Int))
    for i in 1:m
        vecs = Vector{
            PtrArray{
                Int,
                1,
                (1,),
                Tuple{Int},
                Tuple{Nothing},
                Tuple{StrideArrays.StaticInt{1}}
            }
        }(
            undef,
            n
        )
        for j in 1:n
            vecs[j] = alloc(allocator, Vector{Int}, (rand(len),))
            for k in 1:length(vecs[j])
                vecs[j][k] = compute_1(i, j, k)
            end
        end
        s = s + compute_2(vecs)
        clear!(allocator)
    end
    s
end

m, n = 1000, 1000

Random.seed!(0)
do_stuff_1(m, n)

Random.seed!(0)
do_stuff_2(m, n)

Random.seed!(0)
@benchmark do_stuff_1(m, n)

Random.seed!(0)
@benchmark do_stuff_2(m, n)
