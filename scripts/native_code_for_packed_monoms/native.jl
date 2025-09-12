using Revise, InteractiveUtils, Groebner

##############################
## MONOM DIVISIBILITY CHECK ##
##############################

io = open(joinpath(@__DIR__, "native_PackedTuple1_monom_is_divisible.txt"), "w")
code_native(
    io,
    Groebner.monom_is_divisible,
    Tuple{
        Groebner.PackedTuple1{UInt64, UInt8},
        Groebner.PackedTuple1{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

io = open(joinpath(@__DIR__, "native_PackedTuple2_monom_is_divisible.txt"), "w")
code_native(
    io,
    Groebner.monom_is_divisible,
    Tuple{
        Groebner.PackedTuple2{UInt64, UInt8},
        Groebner.PackedTuple2{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

io = open(joinpath(@__DIR__, "native_PackedTuple4_monom_is_divisible.txt"), "w")
code_native(
    io,
    Groebner.monom_is_divisible,
    Tuple{
        Groebner.PackedTuple4{UInt64, UInt8},
        Groebner.PackedTuple4{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

Groebner.monom_is_divisible(
    Groebner.PackedTuple2{UInt64, UInt8}(0x03020100, 0x00000000),
    Groebner.PackedTuple2{UInt64, UInt8}(0x02020100, 0x00000000)
)

# Sasha: probably need to read llvm IR docs to understand why this does not work
# @inline @generated function _packed_vec_uwu(a::NTuple{2,UInt64}, b::NTuple{2,UInt64})
#     N = 16
#     textir = """
#     define i8 @entry(ptr nocapture noundef nonnull readonly align 8 dereferenceable(16) %0, ptr nocapture noundef nonnull readonly align 8 dereferenceable(16) %1) #0 {
#     top:
#         %val.0   = load <2 x i64>, ptr %0
#         %val.1   = load <2 x i64>, ptr %1
#         %av      = bitcast <2 x i64> %val.0 to <$N x i8>
#         %bv      = bitcast <2 x i64> %val.1 to <$N x i8>
#         %mask    = icmp uge <$N x i8> %av, %bv
#         %mask.i  = bitcast <$N x i1> %mask to i$N
#         %res     = icmp eq i$N %mask.i, $(big(2)^N - 1)
#         %retval  = zext i1 %res to i8
#         ret i8 %retval
#     }
#     attributes #0 = { alwaysinline }
#     """
#     quote
#         Base.llvmcall(($textir, "entry"), Bool, Tuple{NTuple{2,UInt64}, NTuple{2,UInt64}}, a, b)
#     end
# end

# _packed_vec_uwu(UInt64.((0x03020100, 0x00000000)), UInt64.((0x02020100, 0x00000000)))
# @code_llvm _packed_vec_uwu(UInt64.((0x03020100, 0x00000000)), UInt64.((0x02020100, 0x00000000)))

##############################
## MONOM CO-PRIMALITY CHECK ##
##############################

io = open(joinpath(@__DIR__, "native_PackedTuple1_monom_is_gcd_const.txt"), "w")
code_native(
    io,
    Groebner.monom_is_gcd_const,
    Tuple{
        Groebner.PackedTuple1{UInt64, UInt8},
        Groebner.PackedTuple1{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

io = open(joinpath(@__DIR__, "native_PackedTuple2_monom_is_gcd_const.txt"), "w")
code_native(
    io,
    Groebner.monom_is_gcd_const,
    Tuple{
        Groebner.PackedTuple2{UInt64, UInt8},
        Groebner.PackedTuple2{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

##############################
##     MONOM LCM COMPUTE    ##
##############################

io = open(joinpath(@__DIR__, "native_PackedTuple1_monom_lcm!.txt"), "w")
code_native(
    io,
    Groebner.monom_lcm!,
    Tuple{
        Groebner.PackedTuple1{UInt64, UInt8},    
        Groebner.PackedTuple1{UInt64, UInt8},
        Groebner.PackedTuple1{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

io = open(joinpath(@__DIR__, "native_PackedTuple2_monom_lcm!.txt"), "w")
code_native(
    io,
    Groebner.monom_lcm!,
    Tuple{
        Groebner.PackedTuple2{UInt64, UInt8},    
        Groebner.PackedTuple2{UInt64, UInt8},
        Groebner.PackedTuple2{UInt64, UInt8}
    },
    debuginfo=:none
)       
close(io)

##############################
##    MONOM HASH COMPUTE    ##
##############################

io = open(joinpath(@__DIR__, "native_PackedTuple1_monom_hash.txt"), "w")
code_native(
    io,
    Groebner.monom_hash,
    Tuple{
        Groebner.PackedTuple1{UInt64, UInt8},    
        Vector{UInt64}
    },
    debuginfo=:none
)       
close(io)

# @inline @generated function _packed_vec_dot2(
#     a::T,
#     b::NTuple{8, MH},
#     ::Type{B}
# ) where {T, MH, B}
#     ts, bs = sizeof(T), sizeof(B)
#     @assert bs * div(ts, bs) == ts
#     epc = div(ts, bs)
#     shift = bs * 8
#     x = :x
#     ans = :(x = $MH(0))
#     for i in 1:epc
#         ans = :($ans;
#         @inbounds iz = $MH(mod(a, $B)) * b[$i];
#         x = x + iz;
#         a = a >> $shift)
#     end
#     :($ans; return $x)
# end

# @inline function _packed_vec_dot3(
#     a::UInt64,
#     b::NTuple{8, UInt64},
# )
#     sum(reinterpret(NTuple{8, UInt8}, a) .* b)
# end
