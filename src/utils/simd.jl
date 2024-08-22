# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# This file provides some fast operations on vectors of integers.

# The use of these functions must be limited.
#
# Debugging LLVM IR is a headache. The LLVM IR textual format is not guaranteed
# to be stable. 
#
# At this moment, we are grudgingly writing the LLVM IR, and we hope that the
# bright future where the compiler does this for us is well within our grasp.

# Functions in this file take their inspiration from
# https://github.com/SciML/FindFirstFunctions.jl

const BitInteger = Union{Int16, Int32, Int64, Int8, UInt16, UInt32, UInt64, UInt8}

function log_simdinfo()
    # avx512f implies avx512bw everywhere except Xeon Phi
    has_avx512 = Bool(has_feature(Val(:x86_64_avx512f)))
    @log :misc """
      CPU (host): $(Symbol(cpu_name()))
      --------------------------
      fma_fast   | $(Bool(fma_fast()))
      registers  | $(Int(register_count()))
      simd width | <$(Int(pick_vector_width(Int32))) x i32>
      has avx512 | $(has_avx512)
      """
    nothing
end

jl_to_llvm_t(::Type{T}) where {T <: BitInteger} = "i$(8*sizeof(T))"

# y, s.t. y >= x and n | y.
align_up(x::Integer, n::Integer) = (x + (n - 1)) & (~(n - 1))

# y, s.t y <= x and n | y.
align_down(x::Integer, n::Integer) = x ⊻ (n - 1)

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

# Permute a part of the array from the given index according to the permutation.
function permute_array!(
    arr::AbstractVector{T},
    perm::Vector{I},
    buf::Vector{T},
    from::Int,
    sz::Int
) where {T, I}
    @invariant length(buf) >= sz && sz <= length(perm)
    @invariant from + sz - 1 <= length(arr)
    @inbounds for i in 1:sz
        buf[i] = arr[perm[i]]
    end
    @inbounds for i in 1:sz
        arr[from + i - 1] = buf[i]
    end
    nothing
end

# Vector functions in this file assume that
# - input vectors are non-negative.
# - input vectors have the same length.

# Returns false if a[i] < b[i] for ANY index i, and true otherwise.
@inline @generated function _vec_not_any_lt(
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
    mask = align_down(typemax(Int), N)
    textir = """
    declare <$N x $llvm_t> @llvm.masked.load.v$(N)$(llvm_t)(<$N x $llvm_t>*, i32, <$N x i1>, <$N x $llvm_t>);
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

# Returns true if a and b are orthogonal, and false otherwise.
@inline @generated function _vec_check_orth(
    a::Vector{T},
    b::Vector{T},
    offset::Int=1
) where {T <: BitInteger}
    N = pick_vector_width_clamp_8(T)

    # Unfortunate case. Default to scalar code.
    if N == 1
        return quote
            @inbounds for j in (1 + offset):length(a)
                if !iszero(a[j]) && !iszero(b[j])
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
    mask = align_down(typemax(Int), N)
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
        %mask1 = icmp ne <$N x $llvm_t> %ai, zeroinitializer
        %mask2 = icmp ne <$N x $llvm_t> %bi, zeroinitializer
        %mask3 = and <$N x i1> %mask1, %mask2
        %compressed = bitcast <$N x i1> %mask3 to i$N
        %matchnotfound1 = icmp eq i$N %compressed, 0
        br i1 %matchnotfound1, label %L30, label %common.ret

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
        %mask11 = icmp ne $llvm_t %savi, 0
        %mask12 = icmp ne $llvm_t %sbvi, 0
        %mask13 = and i1 %mask11, %mask12
        %matchnotfound2 = icmp eq i1 %mask13, 0
        br i1 %matchnotfound2, label %L67, label %common.ret

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

# Returns true if a < b lexicographically, and false otherwise.
@inline @generated function _vec_cmp_lex(
    a::Vector{T},
    b::Vector{T},
    offset::Int=1
) where {T <: BitInteger}
    N = pick_vector_width_clamp_8(T)

    # Unfortunate case. Default to scalar code.
    if N == 1
        return quote
            i = 1 + offset
            @inbounds while i < length(a) && a[i] == b[i]
                i += 1
            end
            @inbounds return a[i] < b[i]
        end
    end

    # The case when IntN exists.
    @assert N in (8, 16, 32, 64)
    B = sizeof(T)
    llvm_t = jl_to_llvm_t(T)
    mask = align_down(typemax(Int), N)
    typecast_iN_to_i64(varname) =
        if N == 64
            ""
        else
            "$(varname)64 = zext i$(N) $(varname)$(N) to i64\n"
        end
    textir = """
    declare i$N @llvm.cttz.i$N(i$N, i1);
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
        %mask = icmp ne <$N x $llvm_t> %ai, %bi
        %compressed = bitcast <$N x i1> %mask to i$N
        %matchfnotound = icmp eq i$N %compressed, 0
        br i1 %matchfnotound, label %L30, label %L17

    L17:
        %tz$N = call i$N @llvm.cttz.i$N(i$N %compressed, i1 true)
        $(typecast_iN_to_i64("%tz"))
        %vis = add nuw i64 %i, %tz64
        %sapi2 = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %vis
        %sbpi2 = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %vis
        %savi2 = load $llvm_t, $llvm_t* %sapi2, align $B
        %sbvi2 = load $llvm_t, $llvm_t* %sbpi2, align $B
        %flag1 = icmp ult $llvm_t %savi2, %sbvi2
        br label %common.ret

    common.ret:
        %retflag = phi i1 [ %flag1, %L17 ], [ 0, %L32 ], [ %flag2, %L51 ], [ 0, %L67 ]
        %retval = zext i1 %retflag to i8
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
        %match = icmp eq $llvm_t %savi, %sbvi
        %flag2 = icmp ult $llvm_t %savi, %sbvi
        br i1 %match, label %L67, label %common.ret

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

# Returns true if a < b reversed lexicographically, and false otherwise.
@inline @generated function _vec_cmp_revlex(
    a::Vector{T},
    b::Vector{T},
    offset::Int=1
) where {T <: BitInteger}
    N = pick_vector_width_clamp_8(T)

    # Unfortunate case. Default to scalar code.
    if N == 1
        return quote
            i = length(a)
            @inbounds while i > 1 + offset && a[i] == b[i]
                i -= 1
            end
            @inbounds return a[i] > b[i]
        end
    end

    # The case when IntN exists.
    @assert N in (8, 16, 32, 64)
    B = sizeof(T)
    llvm_t = jl_to_llvm_t(T)
    mask = align_down(typemax(Int), N)
    typecast_iN_to_i64(varname) =
        if N == 64
            ""
        else
            "$(varname)64 = zext i$(N) $(varname)$(N) to i64\n"
        end

    textir = """
    declare i$N @llvm.ctlz.i$N(i$N, i1);
    define i8 @entry(i64 %0, i64 %1, i64 %2) #0 {
    top:
        %a = inttoptr i64 %0 to $llvm_t*
        %b = inttoptr i64 %1 to $llvm_t*
        %lenm$N = add nsw i64 %2, -$N
        %lenm1 = add nsw i64 %2, -1
        %dosimditer = icmp sge i64 %2, $N
        br i1 %dosimditer, label %L9.lr.ph, label %L32

    L9.lr.ph:
        %len$N = and i64 %2, $mask              ; divisible by N
        %lenmlen$N = sub nsw i64 %lenm1, %len$N
        br label %L9

    L9:
        %i = phi i64 [ %lenm$N, %L9.lr.ph ], [ %vdec, %L30 ]
        %api = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %i
        %bpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %i
        %avi = bitcast $llvm_t* %api to <$N x $llvm_t>*
        %bvi = bitcast $llvm_t* %bpi to <$N x $llvm_t>*
        %ai = load <$N x $llvm_t>, <$N x $llvm_t>* %avi, align $B
        %bi = load <$N x $llvm_t>, <$N x $llvm_t>* %bvi, align $B
        %mask = icmp ne <$N x $llvm_t> %ai, %bi
        %compressed = bitcast <$N x i1> %mask to i$N
        %matchfnotound = icmp eq i$N %compressed, 0
        br i1 %matchfnotound, label %L30, label %L17

    L17:
        %tz$N.rev = call i$N @llvm.ctlz.i$N(i$N %compressed, i1 true)
        %tz$N = sub i$N $(N-1), %tz$N.rev
        $(typecast_iN_to_i64("%tz"))
        %vis = add nsw i64 %i, %tz64
        %sapi2 = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %vis
        %sbpi2 = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %vis
        %savi2 = load $llvm_t, $llvm_t* %sapi2, align $B
        %sbvi2 = load $llvm_t, $llvm_t* %sbpi2, align $B
        %flag1 = icmp ugt $llvm_t %savi2, %sbvi2
        br label %common.ret

    common.ret:
        %retflag = phi i1 [ %flag1, %L17 ], [ 0, %L32 ], [ %flag2, %L51 ], [ 0, %L67 ]
        %retval = zext i1 %retflag to i8
        ret i8 %retval

    L30:
        %vdec = sub nsw i64 %i, $N
        %continuesimd = icmp sge i64 %vdec, 0
        br i1 %continuesimd, label %L9, label %L32

    L32:
        %cumi = phi i64 [ %lenm1, %top ], [ %lenmlen$N, %L30 ]
        %done = icmp eq i64 %cumi, -1
        br i1 %done, label %common.ret, label %L51

    L51:
        %si = phi i64 [ %dec, %L67 ], [ %cumi, %L32 ]
        %sapi = getelementptr inbounds $llvm_t, $llvm_t* %a, i64 %si
        %sbpi = getelementptr inbounds $llvm_t, $llvm_t* %b, i64 %si
        %savi = load $llvm_t, $llvm_t* %sapi, align $B
        %sbvi = load $llvm_t, $llvm_t* %sbpi, align $B
        %match = icmp eq $llvm_t %savi, %sbvi
        %flag2 = icmp ugt $llvm_t %savi, %sbvi
        br i1 %match, label %L67, label %common.ret

    L67:
        %dec = sub nsw i64 %si, 1
        %dobreak = icmp eq i64 %dec, -1
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
