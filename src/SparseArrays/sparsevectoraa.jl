
#=
    The file contains some code from the Julia repository
    https://github.com/JuliaLang/julia

    Namely, we use the source of SparseVectorAA from Julia/SparseArrays in
    our own implementation
=#

import SparseArrays: AbstractSparseVector, nnz, nonzeros, nonzeroinds,
                    findnz, mul!

import Base: size, count, getindex, min, max, map
import Base: +, -, *


struct SparseVectorAA{Tv,Ti<:Integer} <: AbstractSparseVector{Tv,Ti}
    n::Int              # Length of the sparse vector
    nzind::Vector{Ti}   # Indices of stored values
    nzval::Vector{Tv}   # Stored values, typically nonzeros

    zv::Tv              # We add zero value explicitly here

    function SparseVectorAA{Tv,Ti}(n::Integer, nzind::Vector{Ti}, nzval::Vector{Tv}, zv::Tv) where {Tv,Ti<:Integer}
        n >= 0 || throw(ArgumentError("The number of elements must be non-negative."))
        length(nzind) == length(nzval) ||
            throw(ArgumentError("index and value vectors must be the same length"))
        new(convert(Int, n), nzind, nzval, zv)
    end
end


SparseVectorAA(n::Integer, nzind::Vector{Ti}, nzval::Vector{Tv}, zv::Tv) where {Tv,Ti} =
    SparseVectorAA{Tv,Ti}(n, nzind, nzval, zv)


size(x::SparseVectorAA) = (getfield(x, :n), )
count(f, x::SparseVectorAA) = count(f, nonzeros(x)) + f(zero(eltype(x)))*(length(x) - nnz(x))

# implement the nnz - nzrange - nonzeros - rowvals interface for sparse vectors

nnz(x::SparseVectorAA) = length(nonzeros(x))

nonzeros(x::SparseVectorAA) = getfield(x, :nzval)
nonzeroinds(x::SparseVectorAA) = getfield(x, :nzind)

function findnz(x::SparseVectorAA{Tv,Ti}) where {Tv,Ti}
    numnz = nnz(x)

    I = Vector{Ti}(undef, numnz)
    V = Vector{Tv}(undef, numnz)

    nzind = nonzeroinds(x)
    nzval = nonzeros(x)

    @inbounds for i = 1 : numnz
        I[i] = nzind[i]
        V[i] = nzval[i]
    end

    return (I, V)
end

function _spgetindexAA(m::Int, nzind::AbstractVector{Ti}, nzval::AbstractVector{Tv}, i::Integer, zv::Tv) where {Tv,Ti}
    ii = searchsortedfirst(nzind, convert(Ti, i))
    (ii <= m && nzind[ii] == i) ? nzval[ii] : zv
end

function getindex(x::SparseVectorAA, i::Integer)
    checkbounds(x, i)
    _spgetindexAA(nnz(x), nonzeroinds(x), nonzeros(x), i, x.zv)
end


# zero-preserving functions (z->z, nz->nz)
-(x::SparseVectorAA) = SparseVectorAA(length(x), copy(nonzeroinds(x)), -(nonzeros(x)))

# function that does not preserve zeros

macro unarymap_z2nzAA(op, TF)
    esc(quote
        function $(op)(x::AbstractSparseVector{Tv,<:Integer}) where Tv<:$(TF)
            require_one_based_indexing(x)
            v0 = $(op)(x.zv)
            R = typeof(v0)
            xnzind = nonzeroinds(x)
            xnzval = nonzeros(x)
            n = length(x)
            m = length(xnzind)
            y = fill(v0, n)
            @inbounds for j = 1:m
                y[xnzind[j]] = $(op)(xnzval[j])
            end
            y
        end
    end)
end

### Binary Map

# mode:
# 0: f(nz, nz) -> nz, f(z, nz) -> z, f(nz, z) ->  z
# 1: f(nz, nz) -> z/nz, f(z, nz) -> nz, f(nz, z) -> nz
# 2: f(nz, nz) -> z/nz, f(z, nz) -> z/nz, f(nz, z) -> z/nz

function _binarymapAA(f::Function,
                    x::AbstractSparseVector{Tx},
                    y::AbstractSparseVector{Ty},
                    mode::Int) where {Tx,Ty}
    0 <= mode <= 2 || throw(ArgumentError("Incorrect mode $mode."))
    R = typeof(f(x.zv, y.zv))
    n = length(x)
    length(y) == n || throw(DimensionMismatch())

    xnzind = nonzeroinds(x)
    xnzval = nonzeros(x)
    ynzind = nonzeroinds(y)
    ynzval = nonzeros(y)
    mx = length(xnzind)
    my = length(ynzind)
    cap = (mode == 0 ? min(mx, my) : mx + my)::Int

    rind = Vector{Int}(undef, cap)
    rval = Vector{R}(undef, cap)
    ir = 0
    ix = 1
    iy = 1
    zv = x.zv

    ir = (
        mode == 0 ? _binarymap_mode_0AA!(f, mx, my,
            xnzind, xnzval, ynzind, ynzval, rind, rval, zv) :
        mode == 1 ? _binarymap_mode_1AA!(f, mx, my,
            xnzind, xnzval, ynzind, ynzval, rind, rval, zv) :
        _binarymap_mode_2AA!(f, mx, my,
            xnzind, xnzval, ynzind, ynzval, rind, rval, zv)
    )::Int

    resize!(rind, ir)
    resize!(rval, ir)
    return SparseVectorAA(n, rind, rval, zv)
end

function _binarymap_mode_0AA!(f::Function, mx::Int, my::Int,
                            xnzind, xnzval, ynzind, ynzval, rind, rval, zv)
    # f(nz, nz) -> nz, f(z, nz) -> z, f(nz, z) ->  z
    ir = 0; ix = 1; iy = 1
    @inbounds while ix <= mx && iy <= my
        jx = xnzind[ix]
        jy = ynzind[iy]
        if jx == jy
            v = f(xnzval[ix], ynzval[iy])
            ir += 1; rind[ir] = jx; rval[ir] = v
            ix += 1; iy += 1
        elseif jx < jy
            ix += 1
        else
            iy += 1
        end
    end
    return ir
end

function _binarymap_mode_1AA!(f::Function, mx::Int, my::Int,
                            xnzind, xnzval::AbstractVector{Tx},
                            ynzind, ynzval::AbstractVector{Ty},
                            rind, rval, zv) where {Tx,Ty}
    # f(nz, nz) -> z/nz, f(z, nz) -> nz, f(nz, z) -> nz
    ir = 0; ix = 1; iy = 1
    @inbounds while ix <= mx && iy <= my
        jx = xnzind[ix]
        jy = ynzind[iy]
        if jx == jy
            v = f(xnzval[ix], ynzval[iy])
            if v != zv
                ir += 1; rind[ir] = jx; rval[ir] = v
            end
            ix += 1; iy += 1
        elseif jx < jy
            v = f(xnzval[ix], zv)
            ir += 1; rind[ir] = jx; rval[ir] = v
            ix += 1
        else
            v = f(zv, ynzval[iy])
            ir += 1; rind[ir] = jy; rval[ir] = v
            iy += 1
        end
    end
    @inbounds while ix <= mx
        v = f(xnzval[ix], zv)
        ir += 1; rind[ir] = xnzind[ix]; rval[ir] = v
        ix += 1
    end
    @inbounds while iy <= my
        v = f(zv, ynzval[iy])
        ir += 1; rind[ir] = ynzind[iy]; rval[ir] = v
        iy += 1
    end
    return ir
end

function _binarymap_mode_2AA!(f::Function, mx::Int, my::Int,
                            xnzind, xnzval::AbstractVector{Tx},
                            ynzind, ynzval::AbstractVector{Ty},
                            rind, rval, zv) where {Tx,Ty}
    # f(nz, nz) -> z/nz, f(z, nz) -> z/nz, f(nz, z) -> z/nz
    ir = 0; ix = 1; iy = 1
    @inbounds while ix <= mx && iy <= my
        jx = xnzind[ix]
        jy = ynzind[iy]
        if jx == jy
            v = f(xnzval[ix], ynzval[iy])
            if v != zero(v)
                ir += 1; rind[ir] = jx; rval[ir] = v
            end
            ix += 1; iy += 1
        elseif jx < jy
            v = f(xnzval[ix], zv)
            if v != zero(v)
                ir += 1; rind[ir] = jx; rval[ir] = v
            end
            ix += 1
        else
            v = f(zv, ynzval[iy])
            if v != zero(v)
                ir += 1; rind[ir] = jy; rval[ir] = v
            end
            iy += 1
        end
    end
    @inbounds while ix <= mx
        v = f(xnzval[ix], zero(Ty))
        if v != zero(v)
            ir += 1; rind[ir] = xnzind[ix]; rval[ir] = v
        end
        ix += 1
    end
    @inbounds while iy <= my
        v = f(zero(Tx), ynzval[iy])
        if v != zero(v)
            ir += 1; rind[ir] = ynzind[iy]; rval[ir] = v
        end
        iy += 1
    end
    return ir
end

# definition of a few known broadcasted/mapped binary functions — all others defer to HigherOrderFunctions

_bcast_binary_mapAA(f, x, y, mode) = length(x) == length(y) ? _binarymap(f, x, y, mode) : HigherOrderFns._diffshape_broadcast(f, x, y)
for (fun, mode) in [(:+, 1), (:-, 1), (:*, 0), (:min, 2), (:max, 2)]
    fun in (:+, :-) && @eval begin
        # Addition and subtraction can be defined directly on the arrays (without map/broadcast)
        $(fun)(x::SparseVectorAA, y::SparseVectorAA) = _binarymapAA($(fun), x, y, $mode)
    end
    @eval begin
        map(::typeof($fun), x::SparseVectorAA, y::SparseVectorAA) = _binarymapAA($fun, x, y, $mode)
        broadcast(::typeof($fun), x::SparseVectorAA, y::SparseVectorAA) = _bcast_binary_mapAA($fun, x, y, $mode)
    end
end

function mul!(x::SparseVectorAA, y::Re) where {Re<:RingElem}
    map!(c -> c*y, x.nzval, x.nzval)
    x
end

*(x::SparseVectorAA, y) = mul!(deepcopy(x), y)
*(y, x::SparseVectorAA) = mul!(deepcopy(x), y)


function fkeep!(x::SparseVectorAA, f)
    n = length(x::SparseVectorAA)
    nzind = nonzeroinds(x)
    nzval = nonzeros(x)

    x_writepos = 1
    @inbounds for xk in 1:nnz(x)
        xi = nzind[xk]
        xv = nzval[xk]
        # If this element should be kept, rewrite in new position
        if f(xi, xv)
            if x_writepos != xk
                nzind[x_writepos] = xi
                nzval[x_writepos] = xv
            end
            x_writepos += 1
        end
    end

    # Trim x's storage if necessary
    x_nnz = x_writepos - 1
    resize!(nzval, x_nnz)
    resize!(nzind, x_nnz)

    return x
end

dropzeros!(x::SparseVectorAA) = fkeep!(x, (i, x) -> !iszero(x))
