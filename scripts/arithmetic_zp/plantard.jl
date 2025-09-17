using InteractiveUtils

function new_modular(A, B, P::UInt32, R::UInt64)
    C = (((((UInt64(A)*UInt64(B))*R) >>> 32) + 1)*UInt64(P)) >>> 32
    C % UInt32
    # if C == P
    #     0 % UInt32
    # else
    #     C % UInt32
    # end
end

function new_modular_2(A, BR::UInt64, P::UInt32)
    C = ((((UInt64(A)*BR) >>> 32) + 1)*UInt64(P)) >>> 32
    C % UInt32
    # if C == P
    #     0 % UInt32
    # else
    #     C % UInt32
    # end
end

function new_modular_3(C, A, B, P::UInt32, R::UInt64)
    C = (((((C + UInt64(A)*UInt64(B))*R) >>> 32) + 1)*UInt64(P)) >>> 32
    C % UInt32
    # if C == P
    #     0 % UInt32
    # else
    #     C % UInt32
    # end
end

n = 32
phi = (1 + sqrt(5)) / 2
P = UInt32(2^30+3)
@assert P < 2^n / phi
R = UInt64(invmod(P, big(2)^(2n)))

A, B = UInt32(641008720), UInt32(423423)
@assert A < P && B < P

B_repr = UInt64(mod(B*(mod(-mod(BigInt(2)^(2n), P), P)), P))
# B_repr = B

res1 = new_modular(A, B_repr, P, R)

res2 = UInt32((big(A)*B) % P)

res3 = new_modular(res1, 1, P, R)

B_repr_R = B_repr*R

res4 = new_modular_2(A, B_repr_R, P)

@code_native debuginfo=:none new_modular(A, B_repr, P, R)
@code_native debuginfo=:none new_modular_2(A, B_repr_R, P)

@info "" res1 res2 res3 res4

C = UInt64(342342311)
@assert C < P

C_repr = UInt64(mod(C*(mod(-mod(BigInt(2)^(2n), P), P)), P))
res5 = new_modular_3(C_repr, A, B_repr, P, R)

res6 = UInt32((big(C) + big(A)*B) % P)

@info "" res5 res6

# using Groebner, BenchmarkTools, InteractiveUtils
# m = Groebner.SpecializedArithmeticZp(UInt, UInt32, UInt32(2^30+3))
# function uwu(x,y,m)
#     Groebner.mod_p(x*y, m)
# end
# @code_native debuginfo=:none uwu(UInt64(A), UInt64(B_repr), m)
# @btime uwu($(UInt64(A)), $(UInt64(B_repr)), $m)
# @btime new_modular($(UInt64(A)), $(UInt64(B_repr)), $P, $R)
# @btime new_modular_2($(UInt32(A)), $(UInt64(B_repr_R)), $P)

function add_row_1(row, indices, coeffs, val, m)
    @inbounds for i in 1:length(indices)
        ind = indices[i]
        row[ind] = Groebner.mod_p(row[ind] + UInt64(coeffs[i]) * UInt64(val), m)
    end
end

function add_row_2(row, indices, coeffs, val, m)
    @inbounds for i in 1:length(indices)
        ind = indices[i]
        row[ind] = row[ind] + new_modular_2(coeffs[i], val, m)
        row[ind] = row[ind] <= m ? row[ind] : row[ind] - m
    end
end

# @code_native debuginfo=:none add_row_1(zeros(UInt64, 10), [1,3,5], [UInt32(2), UInt32(3), UInt32(4)], UInt32(5), m)
# @code_native debuginfo=:none add_row_2(zeros(UInt64, 10), [1,3,5], [UInt32(2), UInt32(3), UInt32(4)], UInt64(B_repr_R), P)

# @btime add_row_1(__row, __inds, __coeffs, $(UInt32(5)), $m) setup=(__row=rand(UInt64, 1000); __inds=sort(rand(1:1000, 100)); __coeffs=rand(UInt32(1):UInt32(2^30), 100))
# @btime add_row_2(__row, __inds, __coeffs, $(UInt64(B_repr_R)), $P) setup=(__row=rand(UInt64, 1000); __inds=sort(rand(1:1000, 100)); __coeffs=rand(UInt32(1):UInt32(2^30), 100))

#=
.LBB0_4:                                # %L20
        mov     r15, qword ptr [r11 + 8*r8]
        mov     edx, dword ptr [rbx + 4*r8]
        imul    rdx, r10
        add     rdx, qword ptr [rax + 8*r15 - 8]
        mulx    r12, r12, rdi
        shrx    r12, r12, rsi
        imul    r12, rcx
        sub     rdx, r12
        mov     qword ptr [rax + 8*r15 - 8], rdx
        mov     r15, qword ptr [r11 + 8*r8 + 8]
        mov     edx, dword ptr [rbx + 4*r8 + 4]
        imul    rdx, r10
        add     rdx, qword ptr [rax + 8*r15 - 8]
        mulx    r12, r12, rdi
        shrx    r12, r12, rsi
        imul    r12, rcx
        sub     rdx, r12
        mov     qword ptr [rax + 8*r15 - 8], rdx
        add     r8, 2
        cmp     r14, r8
        jne     .LBB0_4

vs.

.LBB0_4:                                # %L20
        mov     r11, qword ptr [rsi + 8*r8]
        mov     ebx, dword ptr [rdi + 4*r8]
        imul    rbx, rcx
        shr     rbx, 32
        inc     rbx
        imul    rbx, rdx
        shr     rbx, 32
        add     qword ptr [rax + 8*r11 - 8], rbx
        mov     r11, qword ptr [rsi + 8*r8 + 8]
        mov     ebx, dword ptr [rdi + 4*r8 + 4]
        imul    rbx, rcx
        shr     rbx, 32
        inc     rbx
        imul    rbx, rdx
        shr     rbx, 32
        add     qword ptr [rax + 8*r11 - 8], rbx
        add     r8, 2
        cmp     r10, r8
        jne     .LBB0_4
=#
