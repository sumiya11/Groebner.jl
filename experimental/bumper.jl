using BenchmarkTools, Bumper

function work0(polys; use_custom_allocator=false)
    if use_custom_allocator
        custom_allocator = Bumper.SlabBuffer()
        @no_escape custom_allocator begin
            work1(polys, custom_allocator)
        end
    else
        work1(polys)
    end
end

# Very important work
function work1(polys, custom_allocator=nothing)
    res = 0
    for poly in polys
        new_poly = work2(poly, custom_allocator)
        res += sum(new_poly)
    end
    res
end

function work2(poly::Vector{T}, ::Nothing) where {T}
    new_poly = Vector{T}(undef, length(poly))
    work3!(new_poly)
end

function work2(poly::Vector{T}, custom_allocator) where {T}
    new_poly = Bumper.alloc!(custom_allocator, T, length(poly))
    work3!(new_poly)
end

function work3!(poly::AbstractVector{T}) where {T}
    poly[1] = one(T)
    for i in 2:length(poly)
        poly[i] = convert(T, i)^3 - poly[i - 1]
    end
    poly
end

###

m, n = 100, 10_000
polys = [rand(UInt32, rand(1:m)) for _ in 1:n];

@btime work0(polys, use_custom_allocator=false)
#   956.800 μs (10001 allocations: 2.55 MiB)
# 0x000000064889383b

@btime work0(polys, use_custom_allocator=true)
#   1.326 ms (5 allocations: 272 bytes)
# 0x000000064889383b
