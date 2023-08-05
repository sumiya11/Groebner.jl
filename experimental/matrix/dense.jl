using SparseArrays, LinearAlgebra
using UnicodePlots
using BenchmarkTools, InteractiveUtils, Profile
using SIMD
using Polyester, CPUSummary
using Base.Threads

macro my_profview(ex)
    :((VSCodeServer.Profile).clear();
    VSCodeServer.Profile.init(n=10^8, delay=0.0000001);
    VSCodeServer.Profile.start_timer();
    $ex;
    VSCodeServer.Profile.stop_timer();
    VSCodeServer.view_profile(;))
end

mutable struct F4matrix{T}
    A::SparseMatrixCSC{T, Int}
    B::SparseMatrixCSC{T, Int}
    C::SparseMatrixCSC{T, Int}
    D::SparseMatrixCSC{T, Int}
    ABCD::SparseMatrixCSC{T, Int}

    uprows::Vector{Vector{Int}}
    lowrows::Vector{Vector{Int}}
    upcoeffs::Vector{Vector{T}}
    lowcoeffs::Vector{Vector{T}}

    pivots::Vector{Vector{Int}}

    result_rows::Vector{Vector{Int}}
    result_coeffs::Vector{Vector{T}}

    upcols::Vector{Vector{Int}}
    upcolscoeffs::Vector{Vector{T}}

    function F4matrix(nup, nlow, nleft, nright; T=UInt64, density=0.01)
        @assert nup == nleft
        @assert nright >= nlow
        m = nup + nlow
        n = nleft + nright
        D_density = 0.3
        C_density = 2 * density
        B_density = 2 * density
        A_density = density / 2
        D = sprand(T, nlow, nright, D_density)
        C = sprand(T, nlow, nleft, C_density)
        B = sprand(T, nup, nright, B_density)
        A = sprand(T, nup, nleft, A_density)
        for i in 2:nup
            A[i, 1:(i - 1)] .= T(0)
        end
        for i in 1:nup
            A[i, i] = T(1)
        end
        dropzeros!(A)
        ABCD = cat(cat(A, B, dims=2), cat(C, D, dims=2), dims=1)
        @debug "Matrix $(m) x $(n) ($(round(100 * nnz(ABCD) / (m*n), digits=2)) % nnz)"

        uprows = Vector{Vector{Int}}(undef, nup)
        lowrows = Vector{Vector{Int}}(undef, nlow)
        upcoeffs = Vector{Vector{T}}(undef, nup)
        lowcoeffs = Vector{Vector{T}}(undef, nlow)
        pivots = Vector{Vector{Int}}(undef, n)

        @inbounds for i in 1:m
            row = ABCD[i, :]
            nzind, nzval = findnz(row)
            if i <= nup
                uprows[i] = nzind
                upcoeffs[i] = nzval
            else
                lowrows[i - nup] = nzind
                lowcoeffs[i - nup] = nzval
            end
        end

        upcols = Vector{Vector{Int}}(undef, nup)
        upcolscoeffs = Vector{Vector{T}}(undef, nup)

        @inbounds for i in 1:nleft
            Ai = @view A[:, i]
            nnz_ind, nnz_val = findnz(Ai)
            upcols[i] = nnz_ind
            upcolscoeffs[i] = nnz_val
        end

        @inbounds for i in 1:nup
            pivots[i] = uprows[i]
        end

        result_rows = Vector{Vector{Int}}(undef, n)
        result_coeffs = Vector{Vector{T}}(undef, n)

        new{T}(
            A,
            B,
            C,
            D,
            ABCD,
            uprows,
            lowrows,
            upcoeffs,
            lowcoeffs,
            pivots,
            result_rows,
            result_coeffs,
            upcols,
            upcolscoeffs
        )
    end
end

function sizes(m::F4matrix)
    size(m.A, 1), size(m.C, 1), size(m.A, 2), size(m.B, 2)
end

function load_sparse_row!(row::Vector{T}, indices, coeffs) where {T}
    @inbounds for i in 1:length(row)
        row[i] = T(0)
    end
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end
    nothing
end

function extract_sparse_row!(indices, coeffs, row::Vector{T}, from::Int, to::Int) where {T}
    z = zero(T)
    j = 1
    @inbounds for i in from:to
        if row[i] != z
            indices[j] = i
            coeffs[j] = row[i]
            j += 1
        end
    end
    nothing
end

function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{Int},
    coeffs::Vector{T}
) where {T}
    @inbounds mul = (typemax(T) - 1) - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j] * coeffs[j]
    end
    nothing
end

function reduce_dense_row_by_dense_row!(row::Vector{T}, b, startcol) where {T}
    @inbounds mul = (typemax(T) - 1) - row[startcol]
    @inbounds for j in startcol:length(row)
        row[j] = row[j] + mul * b[j] * b[j]
    end
    nothing
end

function reduce_low_row!(row, m, i)
    nup, nlow, nleft, nright = sizes(m)
    pivots = m.pivots
    coeffs = m.upcoeffs

    nnz_inds = m.lowrows[i]
    nnz_vals = m.lowcoeffs[i]
    start = nnz_inds[1]
    load_sparse_row!(row, nnz_inds, nnz_vals)

    firstnnz = 0
    nnzcount = 0
    @inbounds for i in start:nleft
        if iszero(row[i])
            continue
        end

        if !isassigned(pivots, i)
            if firstnnz == 0
                firstnnz = i
            end
            nnzcount += 1
            continue
        end

        red_inds = pivots[i]
        red_vals = coeffs[i]
        reduce_dense_row_by_sparse_row!(row, red_inds, red_vals)
    end
    nnzcount == 0, firstnnz, nnzcount
end

function reduce_low_rows_with_upper_rows!(m::F4matrix{T}) where {T}
    nup, nlow, nleft, nright = sizes(m)
    buf = zeros(T, nleft + nright)
    @inbounds for i in 1:nlow
        zeroed, firstnnz, nnzcount = reduce_low_row!(buf, m, i)
        zeroed && continue

        indices = Vector{Int}(undef, nnzcount)
        coeffs = Vector{T}(undef, nnzcount)
        extract_sparse_row!(indices, coeffs, buf, firstnnz, length(buf))

        m.result_rows[i] = indices
        m.result_coeffs[i] = coeffs
    end
    nothing
end

function reduce_low_rows_with_upper_rows_polyester!(m::F4matrix{T}) where {T}
    nup, nlow, nleft, nright = sizes(m)
    @batch threadlocal = zeros(T, nleft + nright) minbatch = 10 for i in 1:nlow
        zeroed, firstnnz, nnzcount = reduce_low_row!(threadlocal, m, i)
        zeroed && continue

        indices = Vector{Int}(undef, nnzcount)
        coeffs = Vector{T}(undef, nnzcount)
        extract_sparse_row!(indices, coeffs, threadlocal, firstnnz, nleft + nright)

        m.result_rows[i] = indices
        m.result_coeffs[i] = coeffs
    end
    nothing
end

function interreduce_upper_rows!(m::F4matrix{T}) where {T}
    nup, nlow, nleft, nright = sizes(m)
    A = m.A
    B = m.B
    row = zeros(T, nright)
    @inbounds for i in nleft:-1:1
        col_nnz_ind = m.upcols[i]
        col_nnz_val = m.upcolscoeffs[i]
        load_sparse_row!(row, nnz_ind, nnz_val)
        for j in 1:(length(nnz_ind) - 1)
            kij = nnz_val[j]
            indices = m.reduce_dense_row_by_sparse_row!(row, indices, coeffs)
            extract_sparse_row!(indices, coeffs, row, 1, 1)
        end
    end
    nothing
end

function benchmark()
    for nlow in (1, 50, 500, 1000, 2000, 5000)
        nright = max(nlow, 100)
        nup = 5000
        nleft = nup

        @info "nleft=$nleft, nright=$nright, nlow=$nlow"

        @info "reduce_low_rows_with_upper_rows!"
        @btime reduce_low_rows_with_upper_rows!(m) setup =
            (m = F4matrix($nup, $nlow, $nleft, $nright, density=0.003))

        @info "reduce_low_rows_with_upper_rows_polyester!"
        @btime reduce_low_rows_with_upper_rows_polyester!(m) setup =
            (m = F4matrix($nup, $nlow, $nleft, $nright, density=0.003))
        # @info "interreduce_upper_rows!"
        # @btime interreduce_upper_rows!(m) setup =
        #     (m = F4matrix($nup, $nlow, $nleft, $nright, density=0.01))
    end
end

benchmark()

m = F4matrix(5000, 1000, 5000, 5000, density=0.003);

# interreduce_upper_rows!(m)

# @my_profview interreduce_upper_rows!(m)

# @code_warntype interreduce_upper_rows!(m)

@time reduce_low_rows_with_upper_rows!(m)

@time reduce_low_rows_with_upper_rows_polyester!(m)

# interreduce_upper_rows!(m)

# @benchmark reduce_low_rows_with_upper_rows!($m)

spy(m.ABCD)
spy(m.A)
