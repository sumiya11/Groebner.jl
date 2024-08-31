# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# This file provides some fast operations on packed vectors of integers.

# checks that a[i] >= b[i] for all i.
@inline @generated function _packed_vec_ge(a::UInt64, b::UInt64)
    N = 8
    textir = """
    define i8 @entry(i64 %0, i64 %1) #0 {
    top:
        %av      = bitcast i64 %0 to <$N x i8>
        %bv      = bitcast i64 %1 to <$N x i8>
        %mask    = icmp uge <$N x i8> %av, %bv
        %mask.i  = bitcast <$N x i1> %mask to i$N
        %res     = icmp eq i$N %mask.i, $(big(2)^N - 1)
        %retval  = zext i1 %res to i8
        ret i8 %retval
    }
    attributes #0 = { alwaysinline }
    """
    quote
        Base.llvmcall(
            ($textir, "entry"),
            Bool,
            Tuple{UInt64, UInt64},
            a,
            b
        )
    end
end

# checks that a, b are orthogonal
@inline @generated function _packed_vec_is_orth(a::UInt64, b::UInt64)
    N = 8
    textir = """
    define i8 @entry(i64 %0, i64 %1) #0 {
    top:
        %av      = bitcast i64 %0 to <$N x i8>
        %bv      = bitcast i64 %1 to <$N x i8>
        %zero    = bitcast i64 0 to <$N x i8>
        %mask.a  = icmp ne <$N x i8> %av, %zero
        %mask.b  = icmp ne <$N x i8> %bv, %zero
        %mask    = and <$N x i1> %mask.a, %mask.b
        %mask.i  = bitcast <$N x i1> %mask to i$N
        %res     = icmp eq i$N %mask.i, 0
        %retval  = zext i1 %res to i8
        ret i8 %retval
    }
    attributes #0 = { alwaysinline }
    """
    quote
        Base.llvmcall(
            ($textir, "entry"),
            Bool,
            Tuple{UInt64, UInt64},
            a,
            b
        )
    end
end

# returns sum a[i]
@inline @generated function _packed_vec_reduce(a::UInt64)
    N = 8
    textir = """
    declare i8 @llvm.vector.reduce.add.v$(N)i8(<$N x i8>)
    define i8 @entry(i64 %0) #0 {
    top:
        %a.v    = bitcast i64 %0 to <$N x i8>
        %sum    = call i8 @llvm.vector.reduce.add.v$(N)i8(<$N x i8> %a.v)
        ret i8 %sum
    }
    attributes #0 = { alwaysinline }
    """
    quote
        Base.llvmcall(
            ($textir, "entry"),
            UInt8,
            Tuple{UInt64},
            a
        )
    end
end

# returns c[i] = max a[i], b[i]
@inline @generated function _packed_vec_max(a::UInt64, b::UInt64)
    N = 8
    textir = """
    declare <$N x i8> @llvm.umax.v$(N)i8(<$N x i8>, <$N x i8>)
    define i64 @entry(i64 %0, i64 %1) #0 {
    top:
        %a.v    = bitcast i64 %0 to <$N x i8>
        %b.v    = bitcast i64 %1 to <$N x i8>
        %max.v  = call <$N x i8> @llvm.umax.v$(N)i8(<$N x i8> %a.v, <$N x i8> %b.v)
        %max    = bitcast <$N x i8> %max.v to i64
        ret i64 %max
    }
    attributes #0 = { alwaysinline }
    """
    quote
        Base.llvmcall(
            ($textir, "entry"),
            UInt64,
            Tuple{UInt64, UInt64},
            a,
            b
        )
    end
end

@inline @generated function _packed_vec_unpack!(
    b::AbstractVector{MH},
    a::T,
    ::Type{B},
    I::Integer
) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    ans = :()
    for i in 1:epc
        ans = :($ans;
        ($i + I > $epc) && return $b;
        @inbounds b[$i] = mod(a, $B);
        a = a >> $shift)
    end
    :($ans; return $b)
end

@inline @generated function _packed_vec_dot(
    a::T,
    b::AbstractVector{MH},
    ::Type{B},
    I::Integer
) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    x = :x
    ans = :(x = $MH(0))
    for i in 1:epc
        ans = :($ans;
        ($i + I > $epc) && return $x;
        @inbounds iz = $MH(mod(a, $B)) * b[$i];
        x = x + iz;
        a = a >> $shift)
    end
    :($ans; return $x)
end
