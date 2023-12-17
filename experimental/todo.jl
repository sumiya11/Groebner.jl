using SIMD, BenchmarkTools

function add_mul!(vres::Vector{Int64}, a::Int32, v::Vector{Int32})
    @inbounds for i in eachindex(vres)
        vres[i] += Int64(a) * Int64(v[i])
    end
end

function add_mul_v2!(vres::Vector{Int64}, a::Int64, v::Vector{Int32})
    @inbounds for i in eachindex(vres)
        vres[i] += a * Int64(v[i])
    end
end
@code_native debuginfo = :none add_mul!([1, 2], Int32(1), Int32[3, 4])

function reduce_dense_row_by_sparse_row_no_remainder_v11!(
    row::Vector{T},
    indices::Vector{I},
    coeffs,
    mul,
    m
) where {I, T}
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] + T(mul) * T(coeffs[j])
        row[idx] = Groebner.mod_p(a, m)
    end

    row
end

function reduce_dense_row_by_sparse_row_no_remainder_v22!(
    row::Vector{T},
    indices::Vector{I},
    coeffs,
    mul,
    m
) where {I, T}
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] + T(mul) * T(coeffs[j])
        row[idx] = a % m
    end

    row
end

function reduce_dense_row_by_sparse_row_no_remainder_2!(
    row::Vector{T},
    indices::Vector{I},
    coeffs,
    mul,
    ::Val{N}
) where {I, T, N}
    @inbounds for j in 1:N:length(indices)
        idx1 = indices[j]
        idx2 = indices[j + 1]
        idx3 = indices[j + 2]
        idx4 = indices[j + 3]
        idx5 = indices[j + 4]
        idx6 = indices[j + 5]
        idx7 = indices[j + 6]
        idx8 = indices[j + 7]

        # c1, c2, c3, c4 = coeffs[j], coeffs[j + 1], coeffs[j + 2], coeffs[j + 3]
        c1, c2 = coeffs[j], coeffs[j + 1]
        c3, c4 = coeffs[j + 2], coeffs[j + 3]
        c5, c6 = coeffs[j + 4], coeffs[j + 5]
        c7, c8 = coeffs[j + 6], coeffs[j + 7]

        # b1, b2, b3, b4 = row[idx1], row[idx2], row[idx3], row[idx4]
        b1, b2 = row[idx1], row[idx2]
        b3, b4 = row[idx3], row[idx4]
        b5, b6 = row[idx5], row[idx6]
        b7, b8 = row[idx7], row[idx8]

        a1 = b1 + T(mul) * T(c1)
        a2 = b2 + T(mul) * T(c2)
        a3 = b3 + T(mul) * T(c3)
        a4 = b4 + T(mul) * T(c4)
        a5 = b5 + T(mul) * T(c5)
        a6 = b6 + T(mul) * T(c6)
        a7 = b7 + T(mul) * T(c7)
        a8 = b8 + T(mul) * T(c8)

        # (a1, a2, a3, a4) = (b1, b2, b3, b4) .+ T(mul) .* (T(c1), T(c2), T(c3), T(c4))
        row[idx1] = a1
        row[idx2] = a2
        row[idx3] = a3
        row[idx4] = a4
        row[idx5] = a5
        row[idx6] = a6
        row[idx7] = a7
        row[idx8] = a8
        # row[idx3] = a3
        # row[idx4] = a4
    end

    row
end

function reduce_dense_row_by_sparse_row_no_remainder_vec!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{C},
    mul::C,
    ::Val{N}
) where {I, T, N, C}
    mul_vec = Vec{N, UInt64}(mul)

    M = div(length(indices), N) * N
    @inbounds for j in 1:N:M
        idx_vec = SIMD.vload(Vec{N, I}, indices, j)

        cfs_vec = SIMD.vload(Vec{N, C}, coeffs, j)
        # cfs_vec_ext = Vec(SIMD.Intrinsics.zext(SIMD.LVec{N, T}, cfs_vec.data))

        row_vec = SIMD.vgather(row, idx_vec)

        # @info "" idx_vec cfs_vec row_vec

        # @info "" cfs_vec_ext

        row_vec_inc = row_vec + mul_vec * cfs_vec

        # @info "" row_vec_inc

        SIMD.vscatter(row_vec_inc, row, idx_vec)
    end

    row
end

function reduce_dense_row_by_sparse_row_no_remainder_vec_2!(
    row::Vector{T},
    indices::Vector{I},
    coeffs::Vector{C},
    mul::C,
    ::Val{N}
) where {I, T, N, C}
    mul_vec = Vec{N, UInt64}(mul)

    M = div(length(indices), 2N) * 2N
    @inbounds for j in 1:N:M
        idx_vec1 = SIMD.vload(Vec{N, I}, indices, j)
        idx_vec2 = SIMD.vload(Vec{N, I}, indices, j)

        cfs_vec1 = SIMD.vload(Vec{N, C}, coeffs, j)
        cfs_vec2 = SIMD.vload(Vec{N, C}, coeffs, j)

        # cfs_vec_ext = Vec(SIMD.Intrinsics.zext(SIMD.LVec{N, T}, cfs_vec.data))

        row_vec1 = SIMD.vgather(row, idx_vec1)
        row_vec2 = SIMD.vgather(row, idx_vec2)

        # @info "" idx_vec cfs_vec row_vec

        # @info "" cfs_vec_ext

        row_vec_inc1 = row_vec1 + mul_vec * cfs_vec1
        row_vec_inc2 = row_vec2 + mul_vec * cfs_vec2

        # @info "" row_vec_inc

        SIMD.vscatter(row_vec_inc1, row, idx_vec1)
        SIMD.vscatter(row_vec_inc2, row, idx_vec2)
    end

    row
end

# @assert reduce_dense_row_by_sparse_row_no_remainder!(
#             UInt64[1, 2, 3, 4],
#             [1, 2, 3, 4],
#             UInt32[3, 1, 0, 10],
#             UInt32(8)
#         ) ==
#         reduce_dense_row_by_sparse_row_no_remainder_vec!(
#             UInt64[1, 2, 3, 4],
#             [1, 2, 3, 4],
#             UInt64[3, 1, 0, 10],
#             UInt64(8),
#             Val(2)
#         ) ==
#         reduce_dense_row_by_sparse_row_no_remainder_2!(
#             UInt64[1, 2, 3, 4],
#             [1, 2, 3, 4],
#             UInt32[3, 1, 0, 10],
#             UInt32(8),
#             Val(2)
#         )

m = Groebner.SpecializedArithmeticZp(UInt64(2^30 + 3))
@code_native debuginfo = :none reduce_dense_row_by_sparse_row_no_remainder_v11!(
    UInt64[1, 2, 3, 4],
    [1, 2, 3, 4],
    UInt64[3, 1, 0, 10],
    UInt64(8),
    m
)

reduce_dense_row_by_sparse_row_no_remainder_vec!(
    UInt64[1, 2, 3, 4, 1],
    [1, 2, 3, 4, 1],
    UInt64[3, 1, 0, 10, 1],
    UInt64(8),
    Val(2)
)

@code_native debuginfo = :none add_mul_v2!(Int64[1, 2, 3, 4], Int64(8), Int32[3, 1, 0, 10])

@code_native debuginfo = :none reduce_dense_row_by_sparse_row_no_remainder_vec_2!(
    UInt64[1, 2, 3, 4],
    [1, 2, 3, 4],
    UInt64[3, 1, 0, 10],
    UInt64(8),
    Val(4)
)

n = 2^10
k = n >> 5
@benchmark reduce_dense_row_by_sparse_row_no_remainder_v11!(v1, i2, v2, c, m) setup = begin
    v1 = rand(UInt, n)
    c = rand(UInt32)
    v2 = rand(UInt32, k)
    i2 = rand(Int(1):Int(n), k)
end

@benchmark reduce_dense_row_by_sparse_row_no_remainder_2!(v1, i2, v2, c, Val(8)) setup =
    begin
        v1 = rand(UInt, n)
        c = rand(UInt32)
        v2 = rand(UInt32, k)
        i2 = rand(Int(1):Int(n), k)
    end

@benchmark reduce_dense_row_by_sparse_row_no_remainder_vec!(v1, i2, v2, c, Val(4)) setup =
    begin
        v1 = rand(UInt, n)
        c = rand(UInt64)
        v2 = rand(UInt64, k)
        i2 = rand(Int(1):Int(n), k)
    end

@benchmark add_mul!(v1, c, v2) setup = begin
    v1 = rand(UInt, n)
    c = rand(UInt32)
    v2 = rand(UInt32, n)
end
