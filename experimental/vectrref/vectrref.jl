using JLD2
using BenchmarkTools
using Plots
using Statistics

name = "matrix_778"
pivs = load((@__DIR__) * "/$name.jld2", "pivs")
matrix = load((@__DIR__) * "/$name.jld2", "matrix")
basis = load((@__DIR__) * "/$name.jld2", "basis")

matrix.nlow, matrix.nup, matrix.ncols
lengths = map(length, pivs[1:(matrix.nup)])
histogram(lengths, bins=100)
sum(lengths), median(lengths), std(lengths), mean(lengths), length(lengths)
quantile(lengths, 0.40),
quantile(lengths, 0.80),
quantile(lengths, 0.90),
quantile(lengths, 0.95)

function rref_1!(heat, row, pivs, coeffs, startcol, n)
    called = []
    @inbounds for i in startcol:n
        iszero(row[i]) && continue
        !isassigned(pivs, i) && continue

        push!(called, i)

        ids = pivs[i]
        cfs = coeffs[i]
        @inbounds for j in 1:length(ids)
            idx = ids[j]
            row[idx] = row[idx] + cfs[j]
            heat[idx] += 1
        end
    end
    called
end

begin
    heat = zeros(Int, matrix.ncols)
    for i in 1:1
        i = 1
        row = zeros(UInt, matrix.ncols)
        Groebner.load_indexed_coefficients!(
            row,
            matrix.lowrows[i],
            basis.coeffs[matrix.low2coef[i]]
        )
        coeffs = Vector{Vector{UInt}}(undef, matrix.ncols)
        for j in 1:length(pivs)
            !isassigned(pivs, j) && continue
            coeffs[j] = basis.coeffs[matrix.up2coef[j]]
        end
        called = rref_1!(heat, row, pivs, coeffs, 1, matrix.ncols)
    end
end
plot(heat)

similarity(a, b) = length(intersect(a, b))
long_ones = filter(x -> length(x) > quantile(lengths, 0.80), pivs[1:(matrix.nup - 1)])
begin
    similarities = []
    for i in 1:(length(long_ones) - 1)
        sim = similarity(long_ones[i], long_ones[i + 1])
        push!(similarities, sim)
    end
    @info "" sum(similarities) sum(length, long_ones)
    plot(similarities)
end
begin
    long_ones2 = sort(long_ones, by=length)
    similarities = []
    for i in 1:(length(long_ones) - 1)
        sim = similarity(long_ones2[i], long_ones[i + 1])
        push!(similarities, sim)
    end
    @info "" sum(similarities) sum(length, long_ones)
    plot(similarities)
end

####

C = UInt64
arithmetic = Groebner.SpecializedBuiltinModularArithmetic(p)
densecoeffs = zeros(C, matrix.ncols)

function reduction(matrix, basis, densecoeffs, magic, pivs)
    @inbounds for i in 1:length(matrix.nlow)
        rowexps = matrix.lowrows[i]
        cfsref = basis.coeffs[matrix.low2coef[i]]
        Groebner.load_indexed_coefficients!(densecoeffs, rowexps, cfsref)
        startcol = rowexps[1]

        zeroed, newrow, newcfs = Groebner.reduce_dense_row_by_known_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            startcol,
            -1,
            magic
        )
        # zeroed = Groebner.reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, Groebner.ColumnIdx(startcol), Groebner.ColumnIdx(-1), magic)
        zeroed && continue

        # matrix.coeffs[i] = newcfs
        # Groebner.normalize_sparse_row!(matrix.coeffs[i], magic)
    end
    nothing
end

@btime for i in 1:8
    reduction($matrix, $basis, $densecoeffs, $magic, $pivs)
end

####

# function reduce_by_pivot!(row::Vector{T}, indices::Vector{Int},
#         cfs::Vector{T}, magic) where {T<:Integer}
#     @inbounds mul = magic.divisor - row[indices[1]]
#     @inbounds for j in 1:length(indices)
#         idx = indices[j]
#         row[idx] = row[idx] + mul * cfs[j]
#     end
#     nothing
# end

# function reduce_by_pivot!(row::Vector{T}, indices::Vector{Int},
#         cfs::Vector{T}, magic) where {T}
#     @inbounds mul = magic.divisor .- row[indices[1]]
#     @inbounds for j in 1:length(indices)
#         idx = indices[j]
#         row[idx] = row[idx] .+ mul .* cfs[j]
#     end
#     nothing
# end

mult_x(arr, X) = [(x, x, x, x) for x in arr]

function setup_1(C, N, K)
    densecoeffs = rand(C, N)
    indices = unique(sort(rand(1:N, K)))
    cfs = rand(C, length(indices))
    densecoeffs, indices, cfs
end

function setup_2(C, N, K)
    densecoeffs = rand(C, N)
    densecoeffs = mult_x(densecoeffs, 4)
    indices = unique(sort(rand(1:N, K)))
    cfs = rand(C, length(indices))
    cfs = mult_x(cfs, 4)
    densecoeffs, indices, cfs
end

(row, indices, cfs) = setup_1(C, N, K)
@code_llvm debuginfo = :none reduce_by_pivot!(row, indices, cfs, magic)

C = UInt64
N = 10000
K = 1000
@btime for i in 1:4
    reduce_by_pivot!(densecoeffs, indices, cfs, $magic)
end setup = ((densecoeffs, indices, cfs) = setup_1(C, N, K))

@btime reduce_by_pivot!(densecoeffs, indices, cfs, $magic) setup =
    ((densecoeffs, indices, cfs) = setup_2(C, N, K))

using BenchmarkTools
using Ghost
using Groebner
using SIMD, VectorizationBase
import Profile

function reduce_by_pivot!(
    row::Vector{T},
    indices::Vector{Int},
    cfs::Vector{T},
    arithmetic
) where {T}
    @inbounds mul = Groebner.divisor(arithmetic) - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = Groebner.mod_x(row[idx] + mul * cfs[j], arithmetic)
    end
    nothing
end

function reduce_by_pivot_ϕ!(
    row::Vector{T},
    indices::Vector{Int},
    cfs::Vector{T},
    arithmetic
) where {T}
    p = Groebner.divisor(arithmetic)
    @inbounds y = p - row[indices[1]]
    ϕ = UInt(div((UInt128(1) << 64) * y, p, Base.RoundUp))
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        x = cfs[j]
        c = Groebner._mul_high(x, ϕ)
        t = row[idx] + x * y - c * p
        row[idx] = ifelse(t > p, t - p, t)
    end
    nothing
end

function reduce_by_pivot_β!(
    row::Vector{T},
    indices::Vector{Int},
    cfs::Vector{T},
    p
) where {T}
    p2 = p^2
    @inbounds c = row[indices[1]]
    m = (length(indices) >> 2) << 2
    @inbounds for j in 1:4:m
        idx0 = indices[j]
        idx1 = indices[j + 1]
        idx2 = indices[j + 2]
        idx3 = indices[j + 3]
        r0, r1, r2, r3 = row[idx0], row[idx1], row[idx2], row[idx3]
        c0, c1, c2, c3 = cfs[j], cfs[j + 1], cfs[j + 2], cfs[j + 3]
        t0 = r0 - c * c0
        t1 = r1 - c * c1
        t2 = r2 - c * c2
        t3 = r3 - c * c3
        t0 = ifelse(!signbit(t0), t0, t0 + p2)
        t1 = ifelse(!signbit(t1), t1, t1 + p2)
        t2 = ifelse(!signbit(t2), t2, t2 + p2)
        t3 = ifelse(!signbit(t3), t3, t3 + p2)
        row[idx0] = t0
        row[idx1] = t1
        row[idx2] = t2
        row[idx3] = t3
    end
    nothing
end

function reduce_by_pivot_β_simd!(
    row::Vector{T},
    indices::Vector{Int},
    cfs::Vector{T},
    p
) where {T}
    p2 = p^2
    N = 8
    @inbounds c = -row[indices[1]]
    mm_p2 = SIMD.Vec{N, T}(p2)
    mm_c = SIMD.Vec{N, T}(c)
    m = length(indices)
    @inbounds for j in 1:N:m
        # mm_idx = SIMD.vload(SIMD.Vec{N, Int}, indices, j)
        mm_y = SIMD.vload(SIMD.Vec{N, T}, cfs, j)
        mm_t = mm_y # SIMD.vgather(row, mm_idx)
        mm_t = mm_t + mm_c * mm_y
        mm_t = mm_t + (mm_t >> UInt8(31)) & mm_p2
        # SIMD.vscatter(mm_t, row, mm_idx)
        SIMD.vstore(mm_t, cfs, j)
    end
    nothing
end

function setup_1(n, k, p, T)
    row = rand(T(0):(T(p) - T(1)), n)
    indices = unique(sort(rand(1:n, k)))
    cfs = rand(T(0):(T(p) - T(1)), length(indices))
    row, indices, cfs
end
function setup_2(n, k, p, T)
    row = rand(T(0):(T(p) - T(1)), n)
    indices = unique(sort(rand(1:n, k)))
    cfs = rand(T(0):(T(p) - T(1)), length(indices))
    row, indices, cfs
end

begin
    p = Int32(32771)
    pu = UInt32(32771)
    arithmetic = Groebner.SpecializedBuiltinModularArithmetic(pu)
    n, k = 2 * 10^1, 2 * 10^1
    r, i, c = setup_2(n, k, p, Int32)
    ru, iu, cu = map(UInt32, r), i, map(UInt32, c)
    r1 = deepcopy(ru)
    r2 = deepcopy(r)
    reduce_by_pivot!(r1, iu, cu, arithmetic)
    reduce_by_pivot_β!(r2, i, c, p)

    map(x -> x % p, r2) .== r1
end

@code_llvm debuginfo = :none reduce_by_pivot!(r1, iu, cu, arithmetic)
@code_llvm debuginfo = :none reduce_by_pivot_β!(r, i, c, p)
@code_llvm debuginfo = :none reduce_by_pivot_β_simd!(r, i, c, p)

@code_native debuginfo = :none reduce_by_pivot_β!(r, i, c, p)

@code_native debuginfo = :none reduce_by_pivot!(r1, iu, cu, arithmetic)
@code_native debuginfo = :none reduce_by_pivot_β!(r, i, c, p)

n, k = 10^3, 10^3
p = UInt32(32771)
arithmetic = Groebner.SpecializedBuiltinModularArithmetic(p)
@btime begin
    reduce_by_pivot!(r, i, c, $arithmetic)
end setup = ((r, i, c) = setup_1($n, $k, $p, UInt32))

p = Int(32771)
@btime begin
    reduce_by_pivot_β!(r, i, c, $p)
end setup = ((r, i, c) = setup_2($n, $k, $p, Int32))

@btime begin
    reduce_by_pivot_β_simd!(r, i, c, $p)
end setup = ((r, i, c) = setup_2($n, $k, $p, Int32))

trace(reduce_by_pivot!, row, indices, cfs, arithmetic)

function rref_1!(row, pivs, coeffs, startcol, n)
    @inbounds for i in startcol:n
        iszero(row[i]) && continue
        !isassigned(pivs, i) && continue

        ids = pivs[i]
        cfs = coeffs[i]
        @inbounds for j in 1:length(ids)
            idx = ids[j]
            row[idx] = row[idx] + cfs[j]
        end
    end
    false
end

function rref_2!(row::Vector{T}, pivs, coeffs, startcol, n, chunksize=3) where {T}
    pivspivs = Vector{Tuple{Int, Vector{Int}, Vector{T}}}(undef, 0)
    buf = zeros(T, chunksize)
    n = n - chunksize
    @inbounds for i in startcol:chunksize:n
        # add pivots
        @inbounds for ii in i:(i + chunksize - 1)
            if !iszero(row[ii])
                # if there is a pivot
                if isassigned(pivs, ii)
                    push!(pivspivs, (1, pivs[ii], coeffs[ii]))
                end
            end
        end
        @inbounds for ii in 1:chunksize

            #
            @inbounds for k in 1:length(pivspivs)
                smth, inds, cfs = pivspivs[k]
                if smth == 0
                    continue
                end
                if inds[smth] > i + ii - 1
                    continue
                end
                buf[ii] += cfs[smth]
                nextsmth = smth < length(inds) ? smth + 1 : 0
                pivspivs[k] = (nextsmth, inds, cfs)
            end

            #
            row[i + ii - 1] += buf[ii]

            #
            buf[ii] = 0
        end
    end
    false
end

function generate_1(n, k, Ti, μ, T)
    pivs, coeffs = Vector{Vector{Int}}(undef, n), Vector{Vector{T}}(undef, n)
    row = Vector{T}(undef, n)
    for i in 1:n
        if rand() < μ
            row[i] = rand(T(1):(T(2)^30))
        else
            row[i] = 0
        end
    end
    idxs = unique(sort(rand(1:n, k)))
    for idx in idxs
        pivs[idx] = unique(sort([idx, rand(idx:n, Ti)...]))
        coeffs[idx] = rand(T(1):(T(2)^30), length(pivs[idx]))
    end
    startcol = 1
    row, pivs, coeffs, startcol, n
end

macro pr(ex, args...)
    return quote
        Profile.clear()
        Profile.init(n=10^7, delay=1e-6)
        Profile.@profile $(esc(ex))
        view_profile(; $(esc.(args)...))
    end
end

begin
    n, k, Ti, μ, T = 1000, 1000, 1000, 0.1, UInt
    row, pivs, coeffs, startcol, n = generate_1(n, k, Ti, μ, T)
    r1, r2 = deepcopy(row), deepcopy(row)

    # @pr for i in 1:100
    #     rref_2!(r1, pivs, coeffs, startcol, n, 3)
    # end

    # @btime rref_2!($r1, $pivs, $coeffs, $startcol, $n, 3)
    @btime rref_1!($r2, $pivs, $coeffs, $startcol, $n)

    count(r1 .== r2)
end

begin
    n = 6
    pivs, coeffs = Vector{Vector{Int}}(undef, n), Vector{Vector{Int}}(undef, n)
    pivs[1] = [1, 2, 4, 5]
    pivs[3] = [3, 4]
    coeffs[1] = [7, 8, 9, 10]
    coeffs[3] = [10, 11]
    row = [1, 1, 1, 1, 1, 1]
    startcol = 1

    rref_1!(row, pivs, coeffs, startcol, n)

    rref_2!(row, pivs, coeffs, startcol, n)

    row
end

using SparseArrays

m, n = 1000, 100000
A = sprand(Float64, m, n, 0.001)
B = sprand(Float64, n, n, 0.001)
B \ A
