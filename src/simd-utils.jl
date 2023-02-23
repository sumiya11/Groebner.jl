# SIMD utils

"""
    always_vectorize()

Always try to vectorize hot loops using SIMD intrinsics.

"""
always_vectorize() = false

"""
    unpack_hi(x::Vec{N, T}, y::Vec{N, T}) where {N, T}

Returns the high part in the product x*y.

"""
@inline function unpack_hi(x::Vec{N, T}, y::Vec{N, T}) where {N, T}
    Vec(_unpack_hi(x.data, y.data))
end

_shuffle_vec(I) = join((string("i32 ", i == :undef ? "undef" : Int32(i::Integer)) for i in I), ", ")

_unpack_hi_bitcast_to(::Type{LVec{2,UInt64}}) = LVec{4,UInt32}
_unpack_hi_bitcast_to(::Type{LVec{4,UInt64}}) = LVec{8,UInt32}
_unpack_hi_shuffle_indices(::Type{LVec{2,UInt64}}) = (1, 4, 3, 5)
_unpack_hi_shuffle_indices(::Type{LVec{4,UInt64}}) = (2, 10, 3, 11, 6, 14, 7, 15)

@generated function _unpack_hi(x::T, y::T) where {T}
    inds = _unpack_hi_shuffle_indices(T)
    shfl = _shuffle_vec(inds)
    U = _unpack_hi_bitcast_to(T)
    M = length(inds)
    llvm = """
        %3 = bitcast $(llvm_type(T)) %0 to $(llvm_type(U))
        %4 = bitcast $(llvm_type(T)) %1 to $(llvm_type(U))
        %5 = shufflevector $(llvm_type(U)) %3, $(llvm_type(U)) %4, <$M x i32> <$shfl>
        %6 = bitcast $(llvm_type(U)) %5 to $(llvm_type(T))
        ret $(llvm_type(T)) %6
        """
    return :(
        $(Expr(:meta, :inline));
        Base.llvmcall($llvm, T, Tuple{T, T}, x, y)
    )
end
