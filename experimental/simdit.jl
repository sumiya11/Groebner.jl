using SIMD
using LoopVectorization
using BenchmarkTools
using CairoMakie: Figure, Axis, Legend, lines!
using CPUSummary
import VectorizationBase
using HostCPUFeatures

function reduce_by_pivot!(row::Vector{T}, indices::Vector{Int},
    cfs::Vector{T}, magic::Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{T}) where {T<:Unsigned}

    @inbounds mul = magic.divisor - row[indices[1]]

    # on our benchmarks usually
    # length(row) / length(indices) varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + mul * cfs[j]) % magic
    end

end

@noinline function do_stuff_noinline(a, p, and, shift)
    a = (a & and) * p
    a = a + a * p
    a << shift
end

@inline function do_stuff_inline(a, p, and, shift)
    a = (a & and) * p
    a = p + a * p
    a << shift
end

function u64_red_inline!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    and = typemax(eltype(row))
    p = magic.divisor
    shift = 7
    @inbounds for i in 1:length(indices)
        idx = indices[i]
        x = row[idx] + mul * cfs[i]
        x = do_stuff_inline(x, p, and, shift)
        row[idx] = x
    end
    row
end

function u64_red_inline_manual!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    N = 2
    and = Vec{N, eltype(row)}(typemax(eltype(row)))
    p = Vec{N, eltype(row)}(magic.divisor)
    shift = Vec{N, Int}(7)
    lane = VecRange{N}(0)
    corr = isodd(length(indices))
    lastidx = length(indices) - corr
    @inbounds for i in 1:N:lastidx
        idx = indices[lane + i]
        x = row[idx] + mul * cfs[lane + i]
        x = do_stuff_inline(x, p, and, shift)
        row[idx] = x
    end
    row
end

function u64_red_inline_manual_4t1!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    N = 4
    T = eltype(row)
    and = Vec{N, T}(typemax(eltype(row)))
    p = Vec{N, T}(magic.divisor)
    shift = Vec{N, Int}(7)
    step = N
    corr = mod(length(indices), step)
    lastidx = length(indices) - corr
    @inbounds for i in 1:step:lastidx
        i1, i2, i3, i4 = i + 0, i + 1, i + 2, i + 3
        idx1, idx2 = indices[i1], indices[i2]
        idx3, idx4 = indices[i3], indices[i4]
        row1234 = Vec{N, T}((row[idx1], row[idx2], row[idx3], row[idx4]))
        cfs1234 = Vec{N, T}((cfs[i1], cfs[i2], cfs[i3], cfs[i4]))
        x1234 = mul*cfs1234
        x1234 = row1234 + x1234
        x1234 = do_stuff_inline(x1234, p, and, shift)
        x1234t = Tuple(x1234)
        row[idx1] = x1234t[1]
        row[idx2] = x1234t[2]
        row[idx3] = x1234t[3]
        row[idx4] = x1234t[4]
    end
    row
end

function u64_red_inline_manual_4t1_2!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    N = 4
    T = eltype(row)
    and = Vec{N, T}(typemax(eltype(row)))
    p = Vec{N, T}(magic.divisor)
    shift = Vec{N, Int}(7)
    step = N
    corr = mod(length(indices), step)
    lastidx = length(indices) - corr
    @inbounds for i in 1:step:lastidx
        # i1, i2, i3, i4 = i + 0, i + 1, i + 2, i + 3
        idx1234 = vload(Vec{N, Int}, indices, i)
        cfs1234 = vload(Vec{N, T}, cfs, i)
        # idx1, idx2 = indices[i1], indices[i2]
        # idx3, idx4 = indices[i3], indices[i4]
        #
        # row1234 = Vec{N, T}((row[idx1], row[idx2], row[idx3], row[idx4]))
        row1234 = vgather(row, idx1234)
        # cfs1234 = Vec{N, T}((cfs[i1], cfs[i2], cfs[i3], cfs[i4]))
        x1234 = mul*cfs1234
        x1234 = row1234 + x1234
        x1234 = do_stuff_inline(x1234, p, and, shift)
        vscatter(x1234, row, idx1234)
        # x1234t = Tuple(x1234)
        # row[idx1] = x1234t[1]
        # row[idx2] = x1234t[2]
        # row[idx3] = x1234t[3]
        # row[idx4] = x1234t[4]
    end
    row
end

function u64_red_inline_manual_4t2!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    N = 4
    unroll = 2
    T = eltype(row)
    step = N*unroll
    and = Vec{N, T}(typemax(eltype(row)))
    p = Vec{N, T}(magic.divisor)
    shift = Vec{N, Int}(7)
    corr = mod(length(indices), step)
    lastidx = length(indices) - corr
    @inbounds for i in 1:step:lastidx
        i1, i2, i3, i4 = i + 0, i + 1, i + 2, i + 3
        i5, i6, i7, i8 = i + 4, i + 5, i + 6, i + 7
        idx1, idx2 = indices[i1], indices[i2]
        idx3, idx4 = indices[i3], indices[i4]
        idx5, idx6 = indices[i5], indices[i6]
        idx7, idx8 = indices[i7], indices[i8]
        row1234 = Vec{N, T}((row[idx1], row[idx2], row[idx3], row[idx4]))
        row5678 = Vec{N, T}((row[idx5], row[idx6], row[idx7], row[idx8]))
        cfs1234 = Vec{N, T}((cfs[i1], cfs[i2], cfs[i3], cfs[i4]))
        cfs5678 = Vec{N, T}((cfs[i5], cfs[i6], cfs[i7], cfs[i8]))
        x1234 = mul*cfs1234
        x5678 = mul*cfs5678
        x1234 = row1234 + x1234
        x5678 = row5678 + x5678
        x1234 = do_stuff_inline(x1234, p, and, shift)
        x5678 = do_stuff_inline(x5678, p, and, shift)
        row[idx1] = x1234[1]
        row[idx2] = x1234[2]
        row[idx3] = x1234[3]
        row[idx4] = x1234[4]
        row[idx5] = x5678[1]
        row[idx6] = x5678[2]
        row[idx7] = x5678[3]
        row[idx8] = x5678[4]
    end
    row
end

function u64_red_inline_manual_8t1!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    N = 8
    T = eltype(row)
    and = Vec{N, T}(typemax(eltype(row)))
    p = Vec{N, T}(magic.divisor)
    shift = Vec{N, Int}(7)
    step = N
    corr = mod(length(indices), step)
    lastidx = length(indices) - corr
    @inbounds for i in 1:step:lastidx
        i1, i2, i3, i4 = i + 0, i + 1, i + 2, i + 3
        i5, i6, i7, i8 = i + 4, i + 5, i + 6, i + 7
        idx1, idx2 = indices[i1], indices[i2]
        idx3, idx4 = indices[i3], indices[i4]
        idx5, idx6 = indices[i5], indices[i6]
        idx7, idx8 = indices[i7], indices[i8]
        idx12345678 = Vec{N, Int}((idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8))
        row12345678 = vgather(row, idx12345678)
        # row12345678 = Vec{N, T}((row[idx1], row[idx2], row[idx3], row[idx4], row[idx5], row[idx6], row[idx7], row[idx8]))
        cfs12345678 = Vec{N, T}((cfs[i1], cfs[i2], cfs[i3], cfs[i4], cfs[i5], cfs[i6], cfs[i7], cfs[i8]))
        x12345678 = mul*cfs12345678
        x12345678 = row12345678 + x12345678
        x12345678 = do_stuff_inline(x12345678, p, and, shift)
        vscatter(x12345678, row, idx12345678)
        # row[idx1] = x12345678[1]
        # row[idx2] = x12345678[2]
        # row[idx3] = x12345678[3]
        # row[idx4] = x12345678[4]
        # row[idx5] = x12345678[5]
        # row[idx6] = x12345678[6]
        # row[idx7] = x12345678[7]
        # row[idx8] = x12345678[8]
    end
    row
end

function u64_red_noinline!(row, indices, cfs, magic)
    @inbounds mul = row[indices[1]]
    and = typemax(eltype(row))
    p = magic.divisor
    shift = 7
    @inbounds for i in eachindex(indices)
        idx = indices[i]
        x = row[idx] + mul * cfs[i]
        x = do_stuff_noinline(x, p, and, shift)
        row[idx] = x
    end
    row
end

function u64_red_baseline!(row, indices, cfs, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + mul * cfs[j]) % magic
    end
    row
end

function u64_red_new!(x::Vector{T}, y::Vector{T}, chi, mu, p, m) where {T}
    @inbounds for j in eachindex(x, y)
        xj = x[j]
        yj = y[j]
        b = (yj * chi) & mu
        c = (yj + widen(b) * p) % T
        x[j] = ifelse(c >= p, c - p, c)
    end
    x
end

function u64_red_new_curious!(x::Vector{T}, y::Vector{T}, chi, mu, p, m) where {T}
    @inbounds for j in eachindex(x, y)
        xj = x[j]
        yj = y[j]
        b = (yj * chi) & mu
        c = yj + widen(b) * p
        d = (c >>> m) % T
        d = d + xj
        x[j] = ifelse(d >= p, d - p, d)
    end
    x
end

function ui_redmod_baseline!(row::Vector{UInt64}, indices, cfs::Vector{UInt64}, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + mul * cfs[j]) % magic
    end
    row
end

function ui_redmod_barrett_fused_half!(row::Vector{UInt64}, indices, cfs::Vector{UInt64}, magic)
    @inbounds mul = row[indices[1]]
    p = magic.divisor
    ϕ = UInt128(div(widen(mul) << (sizeof(UInt64)*8), p, RoundUp))
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = (row[idx] + barrett_mul_mod_half(mul, cfs[j], ϕ, p)) % p
    end
    row
end

function ui_redmod_barrett_fused!(row::Vector{T}, indices, cfs::Vector{T}, magic) where {T}
    @inbounds mul = row[indices[1]]
    p = magic.divisor
    ϕ = (div(widen(mul) << (sizeof(UInt64)*8), p))
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + barrett_mul_mod(mul, cfs[j], ϕ, p)
    end
    row
end

function u64_red_inline_barrett_half!(row, indices, cfs, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    p = magic.divisor
    ϕ = div(mul << 32, p, RoundUp)
    @inbounds for i in 1:length(indices)
        idx = indices[i]
        c = (cfs[i] * ϕ) >> 32
        x = row[idx] + mul*cfs[i] - c*p
        row[idx] = min(x, x - p)
    end
    row
end

function u32_red_inline_barrett_half!(row, indices, cfs, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    p = magic.divisor
    ϕ = div(UInt32(mul) << 32, p, RoundUp)
    @inbounds for i in 1:length(indices)
        idx = indices[i]
        c = ((cfs[i] * ϕ) >> 32) % UInt32
        x = row[idx] + mul*cfs[i] - c*p
        row[idx] = min(x, x - p)
    end
    row
end

function u64_red_inline_barrett_half_turbo!(row, indices, cfs, magic)
    @inbounds mul = magic.divisor - row[indices[1]]
    p = magic.divisor
    ϕ = div(mul << 32, magic)
    @turbo for i in 1:length(indices)
        idx = indices[i]
        c = (cfs[i] * ϕ) >> 32
        x = row[idx] + mul*cfs[i] - c*p
        row[idx] = min(x, x - p)
    end
    row
end

function u64_red_inline_barrett_half_4t32!(row::Vector{T}, indices, cfs, magic) where {T<:UInt32}    
    corr = isodd(length(indices))
    lastidx = length(indices) - corr

    N = 4
    @inbounds mul = magic.divisor - row[indices[1]]
    mulv = SIMD.Vec{N, UInt32}(mul)
    pv = SIMD.Vec{N, UInt32}(magic.divisor)
    ϕv = SIMD.Vec{2, UInt64}(div(UInt64(mul) << 32, magic.divisor))
    zerov = SIMD.Vec{2, UInt64}(0)
    
    #=
    row[idx], cfs[i], mul, ϕ, p
    b, x, y, ϕ, p, s
    
    c = (x * ϕ) >> s
    b + x*y - c*p
    =#
    @inbounds for i in 1:N:lastidx
        cfs1, cfs2, cfs3, cfs4 = cfs[i], cfs[i + 1], cfs[i + 2], cfs[i + 3]
        cfs12 = SIMD.Vec{2, UInt64}((cfs1, cfs2))
        cfs34 = SIMD.Vec{2, UInt64}((cfs3, cfs4))
        
        c12 = cfs12 * ϕv
        c34 = cfs34 * ϕv
        c12 = unpack_hi(c12, zerov)
        c34 = unpack_hi(c34, zerov)

        c1234 = SIMD.Vec{N, UInt32}((c12[1], c12[2], c34[1], c34[2]))

        cfs1234 = SIMD.Vec{4, UInt32}((cfs1, cfs2, cfs3, cfs4))
        idx1234 = vload(SIMD.Vec{N, Int}, indices, i)
        row1234 = vgather(row, idx1234)
        
        x1234 = row1234 + cfs1234 * mulv - c1234*pv
        x1234 = min(x1234, x1234 - pv)
        
        SIMD.vscatter(x1234, row, idx1234)
    end

    row
end

function setup(p, N, T)
    times = Float64[]
    magic = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(
        T(p)
    )
    M = 1000
    times, magic, M
end

function get_times(func, N, T1, T2)
    times, magic, M = setup(2^30+3, N, T1)
    mask = typemax(T2)
    for n in N
        row = rand(T1, M) .& mask
        inds = sort(unique(rand(1:M, n)))
        cfs = rand(T1, length(inds)) .& mask
        bench = @benchmark $func($row, $inds, $cfs, $magic)
        push!(times, minimum(bench.times))
    end
    return times
end

function get_times_prefetch_row(func, N, T)
    times, magic, M = setup(2^31-1, N, T)
    for n in N
        row = rand(T, M)
        inds = sort(unique(rand(1:M, n)))
        cfs = rand(T, length(inds))
        unsafe_prefetch(pointer(row), Val(:write))
        bench = @benchmark $func($row, $inds, $cfs, $magic)
        push!(times, minimum(bench.times))
    end
    return times
end

function get_times_new(func, N, T)
    p = UInt32(2^30 + 3)
    times, magic, M = setup(p, N, T)
    m = 31
    mu = UInt32(2)^m - UInt32(1)

    # Take rho = 1
    rho = UInt32(1)
    # Then,
    chi = div(rho * UInt64(2)^m - 1, p) % UInt32
    for n in N
        row = rand(T, M)
        inds = sort(unique(rand(1:M, n)))
        cfs = rand(T, length(inds))
        bench = @benchmark $func($row, $inds, $cfs, $chi, $mu, $p, $m)
        push!(times, minimum(bench.times))
    end
    return times
end

N = range(32, 256, step=10)

times1 = get_times(ui_redmod_baseline!, N, UInt64, UInt32)
times2 = get_times(u64_red_inline!, N, UInt64, UInt32)
times3 = get_times(u64_red_inline_barrett_half!, N, UInt64, UInt32)

# times2 = get_times(ui_redmod_barrett_fused_half!, N, UInt64)
# times3 = get_times(ui_redmod_barrett_fused!, N, UInt32)
# times2 = get_times(u64_red_inline!, N, UInt64)
# times3 = get_times(u64_red_noinline!, N, UInt64)
# times4 = get_times(u64_red_inline!, N, UInt32)
# times5 = get_times(u64_red_inline_manual_4t1!, N, UInt32)
# times6 = get_times(u64_red_inline_manual_4t1_2!, N, UInt32)
# times6 = get_times(u64_red_inline_manual_8t1!, N, UInt32)
# times7 = get_times(u64_red_inline_manual_4t2!, N, UInt32)
# times8 = get_times_prefetch_row(u64_red_inline_manual_4t1!, N, UInt32)

begin
    fig = Figure(resolution = (1000,500))
    # ax = Axis(fig[1, 1], xlabel="N", ylabel="GUIOPS")
    ax = Axis(fig[1, 1], xlabel="N", ylabel="Runtime (ηs)")

    GUIOPS(N, times) = 5N ./ times

    l1 = lines!(ax, N, times1)
    l2 = lines!(ax, N, times2)
    l3 = lines!(ax, N, times3)
    # l2 = lines!(ax, N, times2)
    # l3 = lines!(ax, N, times3)
    # l4 = lines!(ax, N, times4)
    # l5 = lines!(ax, N, times5)
    # l6 = lines!(ax, N, times6)
    # l7 = lines!(ax, N, times7)
    # l8 = lines!(ax, N, times8)

    fig[1, 2] = Legend(
        fig, [
            l1, 
            l2,
            l3
            # l2, l3, 
            # l4, 
            # l5, 
            # l6, l7,
            # l8
            ], 
        [
            "Baseline (N=UInt64, M=UInt32)", 
            "Something inlined (N=UInt64, M=UInt32)",
            "Barrett half (N=UInt64, M=UInt32)"
            # "Inlined (UInt64)", "@noinline (UInt64)", 
            # "Inlined (UInt32)", 
            # "Inlined Manual 4 × 1 (UInt32)", 
            # "Inlined Manual 8 × 1 (UInt32)",
            # "Inlined Manual 4 × 2 (UInt32)"
            # "Inlined Manual 4 × 1, prefetch row to L1, (UInt32)", 
        ]
    )

    display(fig)
end
