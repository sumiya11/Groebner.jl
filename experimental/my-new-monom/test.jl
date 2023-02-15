# A list of unexhaustive simple tests

include((@__DIR__) * "/my-new-monom.jl")

using Test, TestSetExtensions

# test basics
x = [1, 0, 3, 0, 5]
monom = make_ev(MyNewMonom, x)
@test monom.degrees == UInt64[1, 0, 3, 0, 5]
@test totaldeg(monom) == UInt64(9)

tmp = Vector{UInt64}(undef, length(x))
make_dense!(tmp, monom)
@test tmp == x

# test arithmetic
y = [0, 2, 0, 4, 0]
monom2 = make_ev(MyNewMonom, y)

@test is_gcd_const(monom, monom2)
@test !is_monom_elementwise_eq(monom, monom2)

z = make_zero_ev(MyNewMonom, length(x))
monom_product!(z, monom, monom2)
@test z.degrees == UInt64[1, 2, 3, 4, 5]

@test is_monom_divisible!(z, z, monom2)
@test is_monom_elementwise_eq(z, monom)

# test term order
monom = make_ev(MyNewMonom, x)
monom2 = make_ev(MyNewMonom, y)
@test monom_isless(monom2, monom, Groebner.Lex())

@info "Tests passed"
