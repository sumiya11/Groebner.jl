using Revise, InteractiveUtils, Groebner

##############################
##       BENCHMARK          ##
##############################

using BenchmarkTools

@info "Packed"

@btime Groebner.monom_isless(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000)),
    $(Groebner.DegRevLex())
)
@btime Groebner.monom_product!(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x01010000, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x00000000, 0x00000000, 0x00000000))
)
@btime Groebner.monom_is_divisible(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x03020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000))
)
@btime Groebner.monom_is_gcd_const(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x03020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000))
)
@btime Groebner.monom_lcm!(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x03020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x02020100, 0x00000000, 0x00000000)),
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x00000000, 0x00000000, 0x00000000))
)
@btime Groebner.monom_hash(
    $(Groebner.PackedTuple3{UInt64, UInt8}(0x03020100, 0x00000000, 0x00000000)),
    $(Vector{UInt64}([1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]))
)

@info "Dense"

@btime Groebner.monom_isless(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.DegRevLex())
)
@btime Groebner.monom_product!(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
)
@btime Groebner.monom_is_divisible(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0]))
)
@btime Groebner.monom_is_gcd_const(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0]))
)
@btime Groebner.monom_lcm!(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Groebner.ExponentVector{UInt8}([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
)
@btime Groebner.monom_hash(
    $(Groebner.ExponentVector{UInt8}([2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0])),
    $(Vector{UInt64}([1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]))
)
