using AbstractAlgebra

include("gizmos.jl")

function fukuoka_mq(n::Int, m::Int, char::Int)
    k = GF(char)
    R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    sol = map(_ -> rand(k), 1:n)
    sys = [randpoly(x, 2) for _ in 1:m]
    c = map(f -> evaluate(f, sol), sys)
    sys = sys .- c
    sys
end
