using HostCPUFeatures, InteractiveUtils

const BitInteger = Union{Int16, Int32, Int64, Int8, UInt16, UInt32, UInt64, UInt8}

jl_to_llvm_t(::Type{T}) where {T <: BitInteger} = "i$(8*sizeof(T))"
align_to(x::Integer, N::Integer) = x ⊻ (N - 1)

function llvm_iota(::Type{T}, N::Integer, start::Int=0) where {T <: BitInteger}
    llvm_t = jl_to_llvm_t(T)
    llvm_vec_t = "<$N x $llvm_t>"
    llvm_vec = "<" * join(["$llvm_t $(T(i))" for i in start:(N + start - 1)], ", ") * ">"
    llvm_vec_t, llvm_vec
end

function pick_vector_width_clamp_8(::Type{T}) where {T}
    N = pick_vector_width(T)
    if N in (8, 16, 32)
        return Int(N)
    end
    if N == 64
        return 32
    end
    1
end

# Vector functions in this file assume that
# - input vectors are non-negative.
# - input vectors have the same length.

# Returns false if a[i] < b[i] for ANY index i, and true otherwise.
@inline @generated function _vec_not_any_lt_NEW(
    a::Vector{T},
    b::Vector{T},
    offset::Int=1
) where {T <: BitInteger}
    N = pick_vector_width_clamp_8(T)

    # Unfortunate case. Default to scalar code.
    if N == 1
        return quote
            @inbounds for j in (1 + offset):length(a)
                if a[j] < b[j]
                    return false
                end
            end
            return true
        end
    end

    # The case when IntN exists.
    @assert N in (8, 16, 32, 64)
    B = sizeof(T)
    llvm_t = jl_to_llvm_t(T)
    mask = align_to(typemax(Int), N)
    _, iota = llvm_iota(Int8, N)
    textir = """
    declare <$N x $llvm_t> @llvm.masked.load.v$(N)$(llvm_t)(<$N x $llvm_t>*, i32, <$N x i1>, <$N x $llvm_t>);
    define i8 @entry(i8* %0, i8* %1, i64 %2) #0 {
    top:
        %a = bitcast i8* %0 to $llvm_t*
        %b = bitcast i8* %1 to $llvm_t*
        %lenm$(N-1) = add nsw i64 %2, -$(N-1)
        %dosimditer = icmp ugt i64 %2, $(N-1)
        br i1 %dosimditer, label %L9.lr.ph, label %L32
        L9.lr.ph:
        %len$N = and i64 %2, $mask  ; divisible by N
        br label %L9

    L9:
        %i = phi i64 [ 0, %L9.lr.ph ], [ %vinc, %L30 ]
        %api = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %i
        %bpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %i
        %avi = bitcast $llvm_t* %api to <$N x $llvm_t>*
        %bvi = bitcast $llvm_t* %bpi to <$N x $llvm_t>*
        %ai = load <$N x $llvm_t>, <$N x $llvm_t>* %avi, align $B
        %bi = load <$N x $llvm_t>, <$N x $llvm_t>* %bvi, align $B
        %mask = icmp ult <$N x $llvm_t> %ai, %bi
        %compressed = bitcast <$N x i1> %mask to i$N
        %matchnotfound = icmp eq i$N %compressed, 0
        br i1 %matchnotfound, label %L30, label %common.ret

    common.ret:
        %retval = phi i8 [ 0, %L9 ], [ 1, %L32 ], [ 0, %L51 ], [ 1, %L67 ]
        ret i8 %retval

    L30:
        %vinc = add nuw nsw i64 %i, $N
        %continue = icmp slt i64 %vinc, %lenm$(N-1)
        br i1 %continue, label %L9, label %L32

    L32:
        %cumi = phi i64 [ 0, %top ], [ %len$N, %L30 ]
        %done = icmp eq i64 %cumi, %2
        br i1 %done, label %common.ret, label %L51

    L51:
        %si = phi i64 [ %inc, %L67 ], [ %cumi, %L32 ]
        %sapi = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %si
        %sbpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %si
        %savi = load $llvm_t, $llvm_t* %sapi, align $B
        %sbvi = load $llvm_t, $llvm_t* %sbpi, align $B
        %match = icmp ult $llvm_t %savi, %sbvi
        br i1 %match, label %common.ret, label %L67

    L67:
        %inc = add i64 %si, 1
        %dobreak = icmp eq i64 %inc, %2
        br i1 %dobreak, label %common.ret, label %L51
    }
    attributes #0 = { alwaysinline }
    """
    quote
        GC.@preserve a b begin
            Base.llvmcall(
                ($textir, "entry"),
                Bool,
                Tuple{Ptr{T}, Ptr{T}, Int64},
                pointer(a) + sizeof(T) * offset,
                pointer(b) + sizeof(T) * offset,
                length(a) - offset
            )
        end
    end
end

@inline @generated function _vec_not_any_lt_OLD(
    a::Vector{T},
    b::Vector{T},
    offset::Int=1
) where {T <: BitInteger}
    N = pick_vector_width_clamp_8(T)

    # Unfortunate case. Default to scalar code.
    if N == 1
        return quote
            @inbounds for j in (1 + offset):length(a)
                if a[j] < b[j]
                    return false
                end
            end
            return true
        end
    end

    # The case when IntN exists.
    @assert N in (8, 16, 32, 64)
    B = sizeof(T)
    llvm_t = jl_to_llvm_t(T)
    mask = align_to(typemax(Int), N)
    textir = """
    define i8 @entry(i64 %0, i64 %1, i64 %2) #0 {
    top:
        %a = inttoptr i64 %0 to $llvm_t*
        %b = inttoptr i64 %1 to $llvm_t*
        %lenm$(N-1) = add nsw i64 %2, -$(N-1)
        %dosimditer = icmp ugt i64 %2, $(N-1)
        br i1 %dosimditer, label %L9.lr.ph, label %L32

    L9.lr.ph:
        %len$N = and i64 %2, $mask  ; divisible by N
        br label %L9

    L9:
        %i = phi i64 [ 0, %L9.lr.ph ], [ %vinc, %L30 ]
        %api = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %i
        %bpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %i
        %avi = bitcast $llvm_t* %api to <$N x $llvm_t>*
        %bvi = bitcast $llvm_t* %bpi to <$N x $llvm_t>*
        %ai = load <$N x $llvm_t>, <$N x $llvm_t>* %avi, align $B
        %bi = load <$N x $llvm_t>, <$N x $llvm_t>* %bvi, align $B
        %mask = icmp ult <$N x $llvm_t> %ai, %bi
        %compressed = bitcast <$N x i1> %mask to i$N
        %matchnotfound = icmp eq i$N %compressed, 0
        br i1 %matchnotfound, label %L30, label %common.ret

    common.ret:
        %retval = phi i8 [ 0, %L9 ], [ 1, %L32 ], [ 0, %L51 ], [ 1, %L67 ]
        ret i8 %retval

    L30:
        %vinc = add nuw nsw i64 %i, $N
        %continue = icmp slt i64 %vinc, %lenm$(N-1)
        br i1 %continue, label %L9, label %L32

    L32:
        %cumi = phi i64 [ 0, %top ], [ %len$N, %L30 ]
        %done = icmp eq i64 %cumi, %2
        br i1 %done, label %common.ret, label %L51

    L51:
        %si = phi i64 [ %inc, %L67 ], [ %cumi, %L32 ]
        %sapi = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %si
        %sbpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %si
        %savi = load $llvm_t, $llvm_t* %sapi, align $B
        %sbvi = load $llvm_t, $llvm_t* %sbpi, align $B
        %match = icmp ult $llvm_t %savi, %sbvi
        br i1 %match, label %common.ret, label %L67

    L67:
        %inc = add i64 %si, 1
        %dobreak = icmp eq i64 %inc, %2
        br i1 %dobreak, label %common.ret, label %L51
    }
    attributes #0 = { alwaysinline }
    """
    quote
        GC.@preserve a b begin
            Base.llvmcall(
                ($textir, "entry"),
                Bool,
                Tuple{Ptr{T}, Ptr{T}, Int64},
                pointer(a) + sizeof(T) * offset,
                pointer(b) + sizeof(T) * offset,
                length(a) - offset
            )
        end
    end
end

@inline function foo_old(a::Vector{UInt32})
    textir = """
    define i32 @entry(i64 %0, i64 %1) #0 {
    top:
        %arr = inttoptr i64 %0 to i32*
        %x = getelementptr inbounds i32, i32* %a, i64 0
        ret i32 %x
    }
    attributes #0 = { alwaysinline }
    """
    quote
        GC.@preserve a begin
            Base.llvmcall(($textir, "entry"), Bool, Tuple{Ptr{T}, Int64}, pointer(a), length(a))
        end
    end
end

for i in 1:5
    n = 2^i
    a, b = rand(UInt8, n), rand(UInt8, n)
    res1 = @btime _vec_not_any_lt_OLD($a, $b)
    res2 = @btime _vec_not_any_lt_NEW($a, $b)
    @assert res1 == res2
end

_vec_not_any_lt(UInt32[1, 2, 3, 4, 5, 6, 7, 8, 8], UInt32[1, 2, 3, 4, 5, 6, 7, 8, 9])
@code_llvm _vec_not_any_lt(UInt32[1, 0, 3], UInt32[1, 2, 3])
@code_native _vec_not_any_lt(UInt32[1, 0, 3], UInt32[1, 2, 3])

io = open((@__DIR__) * "/llvm.ll", "w")
code_native(io, _vec_not_any_lt, map(typeof, (UInt16[1, 0, 3], UInt16[1, 2, 3])))
close(io)

#=
L9.lr.ph:
        %len$N = and i64 %2, $mask  ; divisible by N
        br label %L9

    L9:
        %i = phi i64 [ 0, %L9.lr.ph ], [ %vinc, %L30 ]
        %api = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %i
        %bpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %i
        %avi = bitcast $llvm_t* %api to <$N x $llvm_t>*
        %bvi = bitcast $llvm_t* %bpi to <$N x $llvm_t>*
        %ai = load <$N x $llvm_t>, <$N x $llvm_t>* %avi, align $B
        %bi = load <$N x $llvm_t>, <$N x $llvm_t>* %bvi, align $B
        %mask = icmp ult <$N x $llvm_t> %ai, %bi
        %compressed = bitcast <$N x i1> %mask to i$N
        %matchnotfound = icmp eq i$N %compressed, 0
        br i1 %matchnotfound, label %L30, label %common.ret

    common.ret:
        %retval = phi i8 [ 0, %L9 ], [ 1, %L32 ], [ %smatchnotfound.i8, %L51 ]
        ret i8 %retval

    L30:
        %vinc = add nuw nsw i64 %i, $N
        %continue = icmp slt i64 %vinc, %lenm$(N-1)
        br i1 %continue, label %L9, label %L32

    L32:
        %cumi = phi i64 [ 0, %top ], [ %len$N, %L30 ]
        %done = icmp eq i64 %cumi, %2
        br i1 %done, label %common.ret, label %L51

    L51:
        ; %si = phi i64 [ %inc, %L67 ], [ %cumi, %L32 ]
        %lenmod = and i64 %2, $(N-1)  ; mod N
        %lenmod.t = trunc i64 %lenmod to i8
        %lenmod.vec = insertelement <$N x i8> undef, i8 %lenmod.t, i64 0
        %res.si = shufflevector <$N x i8> %lenmod.vec, <$N x i8> undef, <$N x i32> zeroinitializer
        %loadmask = icmp ugt <$N x i8> %res.si, $iota 
        %sapi = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %cumi
        %sbpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %cumi
        %savi = bitcast $llvm_t* %sapi to <$N x $llvm_t>*
        %sbvi = bitcast $llvm_t* %sbpi to <$N x $llvm_t>*
        %sai = call <$N x $llvm_t> @llvm.masked.load.v$(N)$(llvm_t)(<$N x $llvm_t>* %savi, i32 $B, <$N x i1> %loadmask, <$N x $llvm_t> zeroinitializer)
        %sbi = call <$N x $llvm_t> @llvm.masked.load.v$(N)$(llvm_t)(<$N x $llvm_t>* %sbvi, i32 $B, <$N x i1> %loadmask, <$N x $llvm_t> zeroinitializer)
        %smask = icmp ult <$N x $llvm_t> %sai, %sbi
        %scompressed = bitcast <$N x i1> %smask to i$N
        %smatchnotfound = icmp eq i$N %scompressed, 0
        %smatchnotfound.i8 = zext i1 %smatchnotfound to i8
        br label %common.ret
    }
    attributes #0 = { alwaysinline }
=#

########################

@generated function foo_old(a::Vector{UInt32})
    textir = """
    define i32 @entry(i64 %0) #0 {
        %arr = inttoptr i64 %0 to i32*
        %arr.i = getelementptr inbounds i32, i32* %arr, i64 0
        %x = load i32, i32* %arr.i, align 4
        ret i32 %x
    }
    attributes #0 = { alwaysinline }
    """
    quote
        GC.@preserve a begin
            Base.llvmcall(($textir, "entry"), UInt32, Tuple{Ptr{UInt32}}, pointer(a))
        end
    end
end

@generated function foo_new(a::Vector{UInt32})
    textir = """
    define i32 @entry(i8* %0) #0 {
        %arr = bitcast i8* %0 to i32*
        %arr.i = getelementptr inbounds i32, i32* %arr, i64 0
        %x = load i32, i32* %arr.i, align 4
        ret i32 %x
    }
    attributes #0 = { alwaysinline }
    """
    quote
        GC.@preserve a begin
            Base.llvmcall(($textir, "entry"), UInt32, Tuple{Ptr{UInt32}}, pointer(a))
        end
    end
end

@assert foo_old([UInt32(9)]) === foo_new([UInt32(9)]) === UInt32(9)

@code_native foo_old([UInt32(1)])
@code_native foo_new([UInt32(1)])

@code_llvm debuginfo = :none foo_old([UInt32(1)])
@code_llvm debuginfo = :none foo_new([UInt32(1)])

@btime foo_old($([UInt32(9)]));
@btime foo_new($([UInt32(9)]));
