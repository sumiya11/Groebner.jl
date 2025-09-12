using Revise, Groebner, Nemo

R, x = polynomial_ring(QQ, ["x$i" for i in 1:5])
sys = [
    x[1]^2 - 1,
    x[1]^2 - 1,
    x[2]^2 - 2,
    x[2]^2 - 2,
]
@assert dimension(sys) == 3
