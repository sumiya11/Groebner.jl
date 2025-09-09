# https://discourse.julialang.org/t/strange-slowdown-with-threads-greedy-with-bigints/132211

using Base.Threads

function crt(x, y, z)
    M = BigInt(2^3000+3)
    Threads.@threads for i in 1:length(x)
        Base.GMP.MPZ.mul!(x[i], y[i], z[i])
        Base.GMP.MPZ.fdiv_q!(x[i], z[i])
        Base.GMP.MPZ.fdiv_r!(x[i], M)
    end
end

M = 256*256
x = [BigInt() for _ in 1:M]
y = [rand(1:BigInt(2)^3000) for _ in 1:M]
z = [rand(1:BigInt(2)^3000) for _ in 1:M]
@time crt(x, y, z)
@time crt(x, y, z)
