using JLD2
using BenchmarkTools

name = "matrix_9"
pivs = load("experimental/vectrref/$name.jld2", "pivs") 
matrix = load("experimental/vectrref/$name.jld2", "matrix") 
basis = load("experimental/vectrref/$name.jld2", "basis")

####

C = UInt64
magic = Groebner.select_divisor(matrix.coeffs, C(2^31-1))
densecoeffs = zeros(C, matrix.ncols)

function reduction(matrix, basis, densecoeffs, magic, pivs)
    @inbounds for i in 1:length(matrix.nlow)
        rowexps = matrix.lowrows[i]
        cfsref = basis.coeffs[matrix.low2coef[i]]        
        Groebner.load_indexed_coefficients!(densecoeffs, rowexps, cfsref)
        startcol = rowexps[1]
        
        zeroed, newrow, newcfs = Groebner.reduce_dense_row_by_known_pivots_sparse!(densecoeffs, matrix, basis, pivs, startcol, -1, magic)
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

function reduce_by_pivot!(row::Vector{T}, indices::Vector{Int},
        cfs::Vector{T}, magic) where {T<:Integer}
    @inbounds mul = magic.divisor - row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * cfs[j]
    end
    nothing
end

function reduce_by_pivot!(row::Vector{T}, indices::Vector{Int},
        cfs::Vector{T}, magic) where {T}
    @inbounds mul = magic.divisor .- row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] .+ mul .* cfs[j]
    end
    nothing
end

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

(row, indices, cfs)=setup_1(C, N, K)
@code_llvm debuginfo=:none reduce_by_pivot!(row, indices, cfs, magic)

C = UInt64
N = 10000
K = 1000
@btime for i in 1:4
    reduce_by_pivot!(densecoeffs, indices, cfs, $magic)
end setup=((densecoeffs, indices, cfs)=setup_1(C, N, K))

@btime reduce_by_pivot!(densecoeffs, indices, cfs, $magic) setup=((densecoeffs, indices, cfs)=setup_2(C, N, K))

