using AbstractAlgebra
import Nemo

function randpoly(x, d)
    K = base_ring(parent(x[1]))
    M = collect(monomials((sum(x)+1)^d))
    C = map(_ -> rand(K), M)
    f = parent(x[1])(C, exponent_vector.(M, 1))
    f
end

function randsys(n, d)
    R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    [randpoly(x, d) for _ in 1:n]
end

n, d = 3, 5
R, x = polynomial_ring(GF(2^30+3), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
sys = [randpoly(x, d) for _ in 1:n]

n, d = 20, 1
R, x = polynomial_ring(Nemo.GF(2^30+3), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
A = rand(base_ring(parent(x[1])), n, n) 
sys = A*x;

# TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER); @time gb = groebner(sys); show(Groebner._TIMER, allocations=false)
