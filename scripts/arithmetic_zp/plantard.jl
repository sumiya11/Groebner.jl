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

function new_modular_4(C, A, BR, P::UInt32)
    C = ((((C + UInt64(A)*BR) >>> 32) + 1)*UInt64(P)) >>> 32
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

@code_native debuginfo=:none new_modular_3(C_repr, A, B_repr, P, R)

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
    end
end

function add_row_3(row, indices, coeffs, val, m)
    m2 = m^2
    @inbounds for i in 1:length(indices)
        ind = indices[i]
        row[ind] = row[ind] - Int64(coeffs[i]) * Int64(val)
	row[ind] = row[ind] < 0 ? row[ind] + m2 : row[ind] 
    end
end

function add_row_4(row, indices, coeffs, val, m)
    @inbounds for i in 1:length(indices)
        ind = indices[i]
        row[ind] = row[ind] - Int64(coeffs[i]) * Int64(val)
    end
end

# @code_native debuginfo=:none add_row_1(zeros(UInt64, 10), [1,3,5], [UInt32(2), UInt32(3), UInt32(4)], UInt32(5), m)
# @code_native debuginfo=:none add_row_2(zeros(UInt64, 10), [1,3,5], [UInt32(2), UInt32(3), UInt32(4)], UInt64(B_repr_R), P)
@code_native debuginfo=:none add_row_3(zeros(Int64, 10), [1,3,5], [Int32(2), Int32(3), Int32(4)], Int32(1), Int64(P))
@code_native debuginfo=:none add_row_4(zeros(Int64, 10), [1,3,5], [Int32(2), Int32(3), Int32(4)], Int32(1), Int64(P))

using BenchmarkTools, Groebner

M, N = 1000, 1000
m = Groebner.SpecializedArithmeticZp(UInt, UInt32, UInt32(2^30+3)) 
@btime add_row_1(__row, __inds, __coeffs, $(UInt32(5)), $m) setup=(__row=rand(UInt64, M); __inds=sort(rand(1:M, N)); __coeffs=rand(UInt32(1):UInt32(2^30), N))
m = Groebner.SpecializedArithmeticZp(UInt, UInt32, UInt32(2^31-1))
@btime add_row_1(__row, __inds, __coeffs, $(UInt32(5)), $m) setup=(__row=rand(UInt64, M); __inds=sort(rand(1:M, N)); __coeffs=rand(UInt32(1):UInt32(2^30), N))
@btime add_row_2(__row, __inds, __coeffs, $(UInt64(B_repr_R)), $P) setup=(__row=rand(UInt64, M); __inds=sort(rand(1:M, N)); __coeffs=rand(UInt32(1):UInt32(2^30), N))
@btime add_row_3(__row, __inds, __coeffs, $(Int32(5)), $(Int64(P))) setup=(__row=rand(Int64, M); __inds=sort(rand(1:M, N)); __coeffs=rand(Int32(1):Int32(2^30), N))
@btime add_row_4(__row, __inds, __coeffs, $(Int32(5)), $(Int64(P))) setup=(__row=rand(Int64, M); __inds=sort(rand(1:M, N)); __coeffs=rand(Int32(1):Int32(2^30), N))

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

vs.

.LBB0_3:                                # %L21
                                        # =>This Inner Loop Header: Depth=1
        mov     r11, qword ptr [rsi + 8*r10 - 8]
        mov     rbx, qword ptr [rax + 8*r11 - 8]
        movsxd  r14, dword ptr [rdi + 4*r10 - 4]
        imul    r14, rcx
        sub     rbx, r14
        mov     r14, rbx
        sar     r14, 63
        and     r14, r8
        add     r14, rbx
        mov     qword ptr [rax + 8*r11 - 8], r14
        mov     r11, qword ptr [rsi + 8*r10]
        mov     rbx, qword ptr [rax + 8*r11 - 8]
        movsxd  r14, dword ptr [rdi + 4*r10]
        imul    r14, rcx
        sub     rbx, r14
        mov     r14, rbx
        sar     r14, 63
        and     r14, r8
        add     r14, rbx
        mov     qword ptr [rax + 8*r11 - 8], r14
        mov     r11, qword ptr [rsi + 8*r10 + 8]
        mov     rbx, qword ptr [rax + 8*r11 - 8]
        movsxd  r14, dword ptr [rdi + 4*r10 + 4]
        imul    r14, rcx
        sub     rbx, r14
        mov     r14, rbx
        sar     r14, 63
        and     r14, r8
        add     r14, rbx
        mov     qword ptr [rax + 8*r11 - 8], r14
        mov     r11, qword ptr [rsi + 8*r10 + 16]
        mov     rbx, qword ptr [rax + 8*r11 - 8]
        movsxd  r14, dword ptr [rdi + 4*r10 + 8]
        imul    r14, rcx
        sub     rbx, r14
        mov     r14, rbx
        sar     r14, 63
        and     r14, r8
        add     r14, rbx
        mov     qword ptr [rax + 8*r11 - 8], r14
        lea     r11, [r9 + r10]
        add     r11, 4
        add     r10, 4
        cmp     r11, 1
        jne     .LBB0_3
=#

function void()
boot = 100
for _ in 1:boot
        A, B = UInt32(rand(0:P)), UInt32(rand(0:P))
        C = UInt32(rand(0:P))
        B_repr = UInt64(mod(B*(mod(-mod(BigInt(2)^(2n), P), P)), P))
        C_repr = UInt64(mod(C*(mod(-mod(BigInt(2)^(2n), P), P)), P))
        res5 = new_modular_3(C_repr, A, B_repr, P, R)
        res6 = UInt32((big(C) + big(A)*B) % P)
	B_repr_R = B_repr*R
	res7 = new_modular_4(C, A, B_repr_R, P)
	if !(res5 == res6) @info "" res5 res6; @assert false end
	if !(res7 == res6) @info "" res7 res6; @assert false end
end
end
void()

